"""Species-specific, strand-resolved SeiPlant-style GRO-seq models.

The data-preparation path intentionally uses only the Python standard library.
PyTorch, NumPy, SciPy, and Matplotlib are imported only by commands that need
them, so alignment and y-label generation do not depend on the model stack.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
from typing import Iterable, Iterator, Sequence


class ModelWorkflowError(RuntimeError):
    """Raised when inputs or outputs violate a reproducibility contract."""


MANIFEST_FIELDS = (
    "species", "chrom", "start0", "end0", "y_plus", "y_minus", "source"
)


def _open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" \
        else path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_int(value, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _resolve(base: Path, value: str) -> str:
    path = Path(value).expanduser()
    return str((base / path).resolve() if not path.is_absolute() else path.resolve())


def load_model_config(path: str | Path) -> dict:
    """Load, validate, and resolve paths in a strand-model JSON config."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    base = config_path.parent
    config["_config_path"] = str(config_path)
    config["output"] = _resolve(base, config["output"])

    settings = config.setdefault("settings", {})
    settings.setdefault("seed", 42)
    settings.setdefault("window_bp", 1024)
    settings.setdefault("label_step_bp", 512)
    settings.setdefault("negative_ratio", 1.0)
    settings.setdefault("max_ambiguous_fraction", 0.05)
    settings.setdefault("inference_step_bp", 128)
    settings.setdefault("inference_center_bp", 128)
    _require_int(settings["seed"], "settings.seed", 0)
    _require_int(settings["window_bp"], "settings.window_bp")
    _require_int(settings["label_step_bp"], "settings.label_step_bp")
    _require_int(settings["inference_step_bp"], "settings.inference_step_bp")
    _require_int(settings["inference_center_bp"], "settings.inference_center_bp")
    if settings["window_bp"] != 1024:
        raise ValueError("The documented SeiPlant backbone requires 1,024-bp inputs")
    if not isinstance(settings["negative_ratio"], (int, float)) or settings["negative_ratio"] < 0:
        raise ValueError("settings.negative_ratio must be >= 0")
    if not 0 <= settings["max_ambiguous_fraction"] <= 1:
        raise ValueError("settings.max_ambiguous_fraction must be between 0 and 1")
    if settings["inference_center_bp"] > settings["window_bp"]:
        raise ValueError("inference_center_bp cannot exceed window_bp")

    training = config.setdefault("training", {})
    defaults = {
        "batch_size": 256, "epochs": 30, "patience": 5,
        "learning_rate": 1e-5, "weight_decay": 0.0,
        "num_workers": 4, "head_hidden": 256, "amp": True,
    }
    for key, value in defaults.items():
        training.setdefault(key, value)
    for key in ("batch_size", "epochs", "patience", "num_workers", "head_hidden"):
        _require_int(training[key], f"training.{key}", 0 if key == "num_workers" else 1)

    species = config.get("species")
    if not isinstance(species, dict) or not species:
        raise ValueError("species must be a non-empty object")
    slugs: set[str] = set()
    for name, item in species.items():
        if not isinstance(item, dict):
            raise ValueError(f"species.{name} must be an object")
        for key in ("slug", "fasta", "positive_y", "chromosomes"):
            if key not in item:
                raise ValueError(f"species.{name}.{key} is required")
        if item["slug"] in slugs:
            raise ValueError(f"Duplicate species slug: {item['slug']}")
        slugs.add(item["slug"])
        item["fasta"] = _resolve(base, item["fasta"])
        item["positive_y"] = _resolve(base, item["positive_y"])
        if item.get("mask"):
            item["mask"] = _resolve(base, item["mask"])
        splits = item["chromosomes"]
        if set(splits) != {"train", "validation", "test"}:
            raise ValueError(f"{name} chromosomes must define train, validation, and test")
        flat = [str(chrom) for split in splits.values() for chrom in split]
        if len(flat) != len(set(flat)):
            raise ValueError(f"{name} chromosome splits overlap")
        item["chromosomes"] = {key: [str(x) for x in value] for key, value in splits.items()}
    return config


def plan_model(config: dict) -> dict:
    settings = config["settings"]
    return {
        "config": config["_config_path"],
        "output": config["output"],
        "target_columns": ["y_plus", "y_minus"],
        "loss": "mean_squared_error",
        "window_bp": settings["window_bp"],
        "negative_ratio": settings["negative_ratio"],
        "inference": {
            "step_bp": settings["inference_step_bp"],
            "center_bp": settings["inference_center_bp"],
        },
        "species": {
            name: {
                "slug": item["slug"], "fasta": item["fasta"],
                "positive_y": item["positive_y"], "mask": item.get("mask"),
                "chromosomes": item["chromosomes"],
            } for name, item in config["species"].items()
        },
        "training": config["training"],
    }


def read_fai(path: str | Path) -> dict[str, tuple[int, int, int, int]]:
    records: dict[str, tuple[int, int, int, int]] = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ModelWorkflowError(f"Invalid FASTA index row: {line[:120]}")
            records[fields[0]] = tuple(int(x) for x in fields[1:5])
    return records


class IndexedFasta:
    """Small read-only FASTA accessor using the standard samtools .fai format."""

    def __init__(self, fasta: str | Path):
        self.path = Path(fasta)
        self.index = read_fai(str(self.path) + ".fai")
        self._handle = None

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_handle"] = None
        return state

    def _stream(self):
        if self._handle is None:
            self._handle = self.path.open("rb", buffering=1024 * 1024)
        return self._handle

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.index:
            raise KeyError(f"Chromosome {chrom!r} is absent from {self.path}.fai")
        length, offset, line_bases, line_width = self.index[chrom]
        if start < 0 or end < start or end > length:
            raise ValueError(f"Invalid FASTA interval {chrom}:{start}-{end} (length {length})")
        if start == end:
            return ""
        byte_start = offset + (start // line_bases) * line_width + start % line_bases
        last = end - 1
        byte_end = offset + (last // line_bases) * line_width + last % line_bases + 1
        stream = self._stream()
        stream.seek(byte_start)
        return stream.read(byte_end - byte_start).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def read_masks(path: str | None, allowed: set[str]) -> dict[str, list[tuple[int, int]]]:
    masks = {chrom: [] for chrom in allowed}
    if not path:
        return masks
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in allowed:
                continue
            masks[fields[0]].append((int(fields[1]), int(fields[2])))
    for chrom in masks:
        masks[chrom].sort()
    return masks


def overlaps_mask(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    if not intervals:
        return False
    starts = [x[0] for x in intervals]
    index = bisect.bisect_left(starts, end)
    return index > 0 and intervals[index - 1][1] > start


def _ambiguous_fraction(sequence: str) -> float:
    return 1.0 - sum(sequence.count(base) for base in "ATCG") / len(sequence)


def _split_lookup(species: dict) -> dict[str, str]:
    return {chrom: split for split, chroms in species["chromosomes"].items() for chrom in chroms}


def read_positive_labels(path: str | Path, species_name: str, allowed: set[str],
                         window_bp: int) -> list[dict]:
    rows: list[dict] = []
    with _open_text(Path(path)) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"species", "chrom", "start0", "end0", "y_plus", "y_minus"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ModelWorkflowError(f"Missing required y columns in {path}")
        for row in reader:
            if row["chrom"] not in allowed:
                continue
            start, end = int(row["start0"]), int(row["end0"])
            plus, minus = float(row["y_plus"]), float(row["y_minus"])
            if row["species"] != species_name:
                raise ModelWorkflowError(f"Species mismatch in {path}: {row['species']!r}")
            if end - start != window_bp or not (0 <= plus <= 1 and 0 <= minus <= 1):
                raise ModelWorkflowError(f"Invalid y row in {path}: {row}")
            rows.append({"species": species_name, "chrom": row["chrom"],
                         "start0": start, "end0": end, "y_plus": plus,
                         "y_minus": minus, "source": "positive"})
    return rows


def _weighted_chromosome(rng: random.Random, cumulative: list[int], chroms: list[str]) -> str:
    draw = rng.randrange(cumulative[-1])
    return chroms[bisect.bisect_right(cumulative, draw)]


def sample_background(*, rng: random.Random, target: int, chroms: list[str],
                      lengths: dict[str, int], window_bp: int, step_bp: int,
                      occupied: set[tuple[str, int]], masks: dict[str, list[tuple[int, int]]],
                      fasta: IndexedFasta, max_ambiguous_fraction: float,
                      species_name: str) -> tuple[list[dict], dict[str, int]]:
    counts = [max(0, (lengths[chrom] - window_bp) // step_bp + 1) for chrom in chroms]
    usable = [(chrom, count) for chrom, count in zip(chroms, counts) if count]
    if not usable and target:
        raise ModelWorkflowError("No complete windows are available for background sampling")
    chroms = [x[0] for x in usable]
    cumulative: list[int] = []
    total = 0
    for _, count in usable:
        total += count
        cumulative.append(total)
    considered: set[tuple[str, int]] = set()
    rows: list[dict] = []
    rejected = {"positive": 0, "mask": 0, "ambiguous": 0, "duplicate": 0}
    maximum_attempts = max(10000, target * 100)
    attempts = 0
    chrom_order = {chrom: index for index, chrom in enumerate(chroms)}
    while len(rows) < target and attempts < maximum_attempts:
        need = target - len(rows)
        # Draw a reserve, then sort by genomic coordinate before FASTA access.
        # This changes hundreds of thousands of small random disk seeks into a
        # mostly forward scan, which matters for the 14.8-Gb wheat reference.
        batch_target = min(max(need + 10000, math.ceil(need * 1.2)), target * 2 or 1)
        candidates: list[tuple[str, int]] = []
        while len(candidates) < batch_target and attempts < maximum_attempts:
            attempts += 1
            chrom = _weighted_chromosome(rng, cumulative, chroms)
            n_windows = (lengths[chrom] - window_bp) // step_bp + 1
            start = rng.randrange(n_windows) * step_bp
            key = (chrom, start)
            if key in occupied:
                rejected["positive"] += 1
                continue
            if key in considered:
                rejected["duplicate"] += 1
                continue
            considered.add(key)
            if overlaps_mask(masks[chrom], start, start + window_bp):
                rejected["mask"] += 1
                continue
            candidates.append(key)
        candidates.sort(key=lambda key: (chrom_order[key[0]], key[1]))
        for chrom, start in candidates:
            sequence = fasta.fetch(chrom, start, start + window_bp)
            if len(sequence) != window_bp or _ambiguous_fraction(sequence) > max_ambiguous_fraction:
                rejected["ambiguous"] += 1
                continue
            rows.append({"species": species_name, "chrom": chrom, "start0": start,
                         "end0": start + window_bp, "y_plus": 0.0,
                         "y_minus": 0.0, "source": "background"})
            if len(rows) == target:
                break
    if len(rows) != target:
        raise ModelWorkflowError(
            f"Only sampled {len(rows)}/{target} background windows after {attempts} attempts"
        )
    rejected["attempts"] = attempts
    return rows, rejected


def _write_manifest(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(path, "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sequence_codes(sequence: str) -> bytes:
    table = getattr(_sequence_codes, "_table", None)
    if table is None:
        values = bytearray([4] * 256)
        for base, code in ((b"A", 0), (b"T", 1), (b"C", 2), (b"G", 3)):
            values[base[0]] = code
        table = bytes(values)
        _sequence_codes._table = table
    return sequence.upper().encode("ascii").translate(table)


def _write_sequence_cache(path: Path, rows: list[dict], fasta: IndexedFasta,
                          window_bp: int) -> None:
    """Write one byte/base in manifest order for fast memory-mapped training."""
    temporary = path.with_name(path.name + ".building")
    with temporary.open("wb", buffering=4 * 1024 * 1024) as stream:
        for row in rows:
            sequence = fasta.fetch(row["chrom"], row["start0"], row["end0"])
            if len(sequence) != window_bp:
                raise ModelWorkflowError(f"Short FASTA sequence for {row}")
            stream.write(_sequence_codes(sequence))
    expected = len(rows) * window_bp
    if temporary.stat().st_size != expected:
        raise ModelWorkflowError(f"Sequence-cache size mismatch: {temporary}")
    temporary.replace(path)


def select_species(config: dict, requested: Sequence[str] | None) -> list[tuple[str, dict]]:
    if not requested:
        return list(config["species"].items())
    wanted = set(requested)
    selected = [(name, item) for name, item in config["species"].items()
                if name in wanted or item["slug"] in wanted]
    found = {name for name, _ in selected} | {item["slug"] for _, item in selected}
    missing = wanted - found
    if missing:
        raise ValueError(f"Unknown species: {', '.join(sorted(missing))}")
    return selected


def prepare_manifests(config: dict, requested: Sequence[str] | None = None,
                      force: bool = False) -> dict:
    """Build chromosome-separated positive/background manifests."""
    output = Path(config["output"])
    settings = config["settings"]
    summaries = {}
    for species_index, (name, item) in enumerate(select_species(config, requested)):
        destination = output / item["slug"] / "manifests"
        summary_path = destination / "summary.json"
        if summary_path.exists() and not force:
            raise ModelWorkflowError(f"Refusing to overwrite {summary_path}; pass --force")
        fasta_path = Path(item["fasta"])
        for required in (fasta_path, Path(str(fasta_path) + ".fai"), Path(item["positive_y"])):
            if not required.is_file():
                raise ModelWorkflowError(f"Missing input: {required}")
        fasta = IndexedFasta(fasta_path)
        try:
            lengths = {chrom: values[0] for chrom, values in fasta.index.items()}
            split_for = _split_lookup(item)
            missing = set(split_for) - set(lengths)
            if missing:
                raise ModelWorkflowError(f"Chromosomes absent from FASTA for {name}: {sorted(missing)}")
            masks = read_masks(item.get("mask"), set(split_for))
            positives = read_positive_labels(item["positive_y"], name, set(split_for), settings["window_bp"])
            kept: dict[str, list[dict]] = {x: [] for x in ("train", "validation", "test")}
            dropped = {"mask": 0, "ambiguous": 0}
            occupied: dict[str, set[tuple[str, int]]] = {x: set() for x in kept}
            for row in positives:
                split = split_for[row["chrom"]]
                start, end = row["start0"], row["end0"]
                if overlaps_mask(masks[row["chrom"]], start, end):
                    dropped["mask"] += 1
                    continue
                sequence = fasta.fetch(row["chrom"], start, end)
                if len(sequence) != settings["window_bp"] or \
                        _ambiguous_fraction(sequence) > settings["max_ambiguous_fraction"]:
                    dropped["ambiguous"] += 1
                    continue
                kept[split].append(row)
                occupied[split].add((row["chrom"], start))

            split_summary = {}
            seed = settings["seed"] + species_index * 1000
            for split_index, split in enumerate(("train", "validation", "test")):
                target = int(round(len(kept[split]) * settings["negative_ratio"]))
                background, rejected = sample_background(
                    rng=random.Random(seed + split_index), target=target,
                    chroms=item["chromosomes"][split], lengths=lengths,
                    window_bp=settings["window_bp"], step_bp=settings["label_step_bp"],
                    occupied=occupied[split], masks=masks, fasta=fasta,
                    max_ambiguous_fraction=settings["max_ambiguous_fraction"],
                    species_name=name,
                )
                rows = kept[split] + background
                chromosome_order = {chrom: index for index, chrom in enumerate(item["chromosomes"][split])}
                rows.sort(key=lambda row: (chromosome_order[row["chrom"]], row["start0"], row["source"]))
                manifest = destination / f"{split}.tsv.gz"
                _write_manifest(manifest, rows)
                sequence_cache = destination / f"{split}.sequences.uint8"
                _write_sequence_cache(sequence_cache, rows, fasta, settings["window_bp"])
                split_summary[split] = {
                    "positive": len(kept[split]), "background": len(background),
                    "total": len(rows), "background_rejections": rejected,
                    "manifest": str(manifest), "sha256": _sha256(manifest),
                    "sequence_cache": str(sequence_cache),
                    "sequence_cache_sha256": _sha256(sequence_cache),
                    "sequence_cache_encoding": "one unsigned byte/base: A=0,T=1,C=2,G=3,other=4",
                }
            summary = {
                "schema_version": 1, "created_unix": int(time.time()),
                "species": name, "slug": item["slug"], "fasta": str(fasta_path),
                "fasta_fai_sha256": _sha256(Path(str(fasta_path) + ".fai")),
                "positive_y": item["positive_y"],
                "positive_y_sha256": _sha256(Path(item["positive_y"])),
                "mask": item.get("mask"), "settings": settings,
                "chromosomes": item["chromosomes"], "dropped_positive": dropped,
                "splits": split_summary,
            }
            destination.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            summaries[name] = summary
        finally:
            fasta.close()
    return summaries


def _lazy_model_imports():
    try:
        import numpy as np
        import torch
        import torch.nn as nn
        from scipy.interpolate import BSpline
    except ImportError as exc:
        raise ModelWorkflowError(
            "Model dependencies are missing; create and activate environment-model.yml"
        ) from exc
    return np, torch, nn, BSpline


def spline_basis(length: int, degrees_of_freedom: int):
    np, _, _, BSpline = _lazy_model_imports()
    degree = 3
    n_inner = degrees_of_freedom - (degree + 1)
    inner = np.linspace(0, length - 1, n_inner + 2)[1:-1] if n_inner > 0 else np.array([])
    knots = np.concatenate((np.repeat(0.0, degree + 1), inner,
                            np.repeat(float(length - 1), degree + 1)))
    x = np.arange(length, dtype=float)
    basis = np.empty((length, degrees_of_freedom), dtype=np.float32)
    for index in range(degrees_of_freedom):
        coefficients = np.zeros(degrees_of_freedom)
        coefficients[index] = 1.0
        basis[:, index] = BSpline(knots, coefficients, degree)(x)
    return basis


def build_model(head_hidden: int = 256):
    """Construct the Figure-1A-style 1,024-bp, two-output network."""
    _, torch, nn, _ = _lazy_model_imports()

    class BSplineTransformation(nn.Module):
        def __init__(self, df: int = 16):
            super().__init__()
            self.df = df
            self.register_buffer("basis", torch.empty(0), persistent=True)

        def forward(self, x):
            if self.basis.numel() == 0 or self.basis.shape[0] != x.shape[-1]:
                self.basis = torch.from_numpy(spline_basis(x.shape[-1], self.df)).to(x.device)
            return torch.matmul(x, self.basis)

    class StrandSeiPlant(nn.Module):
        def __init__(self):
            super().__init__()
            self.lconv1 = nn.Sequential(nn.Conv1d(4, 480, 9, padding=4),
                                        nn.Conv1d(480, 480, 9, padding=4))
            self.conv1 = nn.Sequential(nn.Conv1d(480, 480, 9, padding=4), nn.ReLU(inplace=True),
                                       nn.Conv1d(480, 480, 9, padding=4), nn.ReLU(inplace=True))
            self.lconv2 = nn.Sequential(nn.MaxPool1d(4, 4), nn.Dropout(0.2),
                                        nn.Conv1d(480, 640, 9, padding=4),
                                        nn.Conv1d(640, 640, 9, padding=4))
            self.conv2 = nn.Sequential(nn.Dropout(0.2), nn.Conv1d(640, 640, 9, padding=4),
                                       nn.ReLU(inplace=True), nn.Conv1d(640, 640, 9, padding=4),
                                       nn.ReLU(inplace=True))
            self.lconv3 = nn.Sequential(nn.MaxPool1d(4, 4), nn.Dropout(0.2),
                                        nn.Conv1d(640, 960, 9, padding=4),
                                        nn.Conv1d(960, 960, 9, padding=4))
            self.conv3 = nn.Sequential(nn.Dropout(0.2), nn.Conv1d(960, 960, 9, padding=4),
                                       nn.ReLU(inplace=True), nn.Conv1d(960, 960, 9, padding=4),
                                       nn.ReLU(inplace=True))
            self.dilated = nn.ModuleList([
                nn.Sequential(nn.Dropout(0.1), nn.Conv1d(960, 960, 5, dilation=d, padding=2*d),
                              nn.ReLU(inplace=True)) for d in (2, 4, 8, 16, 25)
            ])
            self.spline = BSplineTransformation(16)
            self.head = nn.Sequential(nn.Dropout(0.5), nn.Linear(960 * 16, head_hidden),
                                      nn.ReLU(inplace=True), nn.Linear(head_hidden, 2), nn.Sigmoid())

        def forward(self, x):
            l1 = self.lconv1(x); x = self.conv1(l1) + l1
            l2 = self.lconv2(x); x = self.conv2(l2) + l2
            l3 = self.lconv3(x); x = self.conv3(l3) + l3
            for block in self.dilated:
                x = x + block(x)
            x = self.spline(x).flatten(1)
            return self.head(x)

    return StrandSeiPlant()


def one_hot(sequence: str):
    np, _, _, _ = _lazy_model_imports()
    lookup = getattr(one_hot, "_lookup", None)
    if lookup is None:
        lookup = np.full(256, -1, dtype=np.int8)
        for base, channel in ((b"A", 0), (b"T", 1), (b"C", 2), (b"G", 3)):
            lookup[base[0]] = channel
        one_hot._lookup = lookup
    encoded = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
    channels = lookup[encoded]
    valid = channels >= 0
    positions = np.nonzero(valid)[0]
    array = np.zeros((4, len(sequence)), dtype=np.float32)
    array[channels[valid], positions] = 1.0
    return array


def one_hot_codes(codes):
    np, _, _, _ = _lazy_model_imports()
    codes = np.asarray(codes, dtype=np.uint8)
    valid = codes < 4
    positions = np.nonzero(valid)[0]
    array = np.zeros((4, len(codes)), dtype=np.float32)
    array[codes[valid], positions] = 1.0
    return array


def _read_manifest(path: Path) -> list[dict]:
    with _open_text(path) as stream:
        rows = []
        for row in csv.DictReader(stream, delimiter="\t"):
            rows.append({"chrom": row["chrom"], "start0": int(row["start0"]),
                         "end0": int(row["end0"]), "y_plus": float(row["y_plus"]),
                         "y_minus": float(row["y_minus"])})
    return rows


def _dataset_class():
    _, torch, _, _ = _lazy_model_imports()

    class GenomeDataset(torch.utils.data.Dataset):
        def __init__(self, manifest: Path, sequence_cache: Path):
            self.rows = _read_manifest(manifest)
            self.sequence_cache = sequence_cache
            if not self.rows:
                raise ModelWorkflowError(f"Empty training manifest: {manifest}")
            self.window_bp = self.rows[0]["end0"] - self.rows[0]["start0"]
            expected = len(self.rows) * self.window_bp
            if not sequence_cache.is_file() or sequence_cache.stat().st_size != expected:
                raise ModelWorkflowError(f"Missing or invalid sequence cache: {sequence_cache}")
            self._cache = None

        def __getstate__(self):
            state = self.__dict__.copy(); state["_cache"] = None
            return state

        def _sequences(self):
            if self._cache is None:
                import numpy as np
                self._cache = np.memmap(self.sequence_cache, dtype=np.uint8, mode="r",
                                        shape=(len(self.rows), self.window_bp))
            return self._cache

        def __len__(self):
            return len(self.rows)

        def __getitem__(self, index):
            row = self.rows[index]
            x = torch.from_numpy(one_hot_codes(self._sequences()[index]))
            y = torch.tensor([row["y_plus"], row["y_minus"]], dtype=torch.float32)
            return x, y

    return GenomeDataset


def _seed_everything(seed: int) -> None:
    np, torch, _, _ = _lazy_model_imports()
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def _metrics(observed, predicted) -> dict:
    np, _, _, _ = _lazy_model_imports()
    from scipy.stats import spearmanr
    result = {}
    for index, strand in enumerate(("plus", "minus")):
        y, p = observed[:, index], predicted[:, index]
        pearson = float(np.corrcoef(y, p)[0, 1]) if np.std(y) and np.std(p) else None
        spearman = float(spearmanr(y, p).statistic) if np.unique(y).size > 1 and np.unique(p).size > 1 else None
        result[strand] = {
            "n": int(len(y)), "mse": float(np.mean((y - p) ** 2)),
            "mae": float(np.mean(np.abs(y - p))), "pearson_r": pearson,
            "spearman_rho": spearman,
        }
    return result


def _evaluate(model, loader, device):
    np, torch, _, _ = _lazy_model_imports()
    model.eval(); observed = []; predicted = []
    with torch.no_grad():
        for x, y in loader:
            prediction = model(x.to(device, non_blocking=True)).cpu().numpy()
            predicted.append(prediction); observed.append(y.numpy())
    return np.concatenate(observed), np.concatenate(predicted)


def train_species(config: dict, requested: Sequence[str] | None = None,
                  device_name: str | None = None, force: bool = False) -> dict:
    np, torch, nn, _ = _lazy_model_imports()
    settings, training = config["settings"], config["training"]
    results = {}
    for name, item in select_species(config, requested):
        _seed_everything(settings["seed"])
        root = Path(config["output"]) / item["slug"]
        model_dir = root / "model"
        best_path = model_dir / "best.pt"
        if best_path.exists() and not force:
            raise ModelWorkflowError(f"Refusing to overwrite {best_path}; pass --force")
        manifests = root / "manifests"
        for split in ("train", "validation", "test"):
            if not (manifests / f"{split}.tsv.gz").is_file():
                raise ModelWorkflowError(f"Missing manifest; run prepare first: {split}")
        device = torch.device(device_name or ("cuda" if torch.cuda.is_available() else "cpu"))
        GenomeDataset = _dataset_class()
        datasets = {split: GenomeDataset(manifests / f"{split}.tsv.gz",
                                         manifests / f"{split}.sequences.uint8")
                    for split in ("train", "validation", "test")}
        generator = torch.Generator().manual_seed(settings["seed"])
        loaders = {
            "train": torch.utils.data.DataLoader(
                datasets["train"], batch_size=training["batch_size"], shuffle=True,
                num_workers=training["num_workers"], pin_memory=device.type == "cuda",
                persistent_workers=training["num_workers"] > 0, generator=generator),
            "validation": torch.utils.data.DataLoader(
                datasets["validation"], batch_size=training["batch_size"], shuffle=False,
                num_workers=training["num_workers"], pin_memory=device.type == "cuda"),
            "test": torch.utils.data.DataLoader(
                datasets["test"], batch_size=training["batch_size"], shuffle=False,
                num_workers=training["num_workers"], pin_memory=device.type == "cuda"),
        }
        model = build_model(training["head_hidden"]).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=training["learning_rate"],
                                     weight_decay=training["weight_decay"])
        loss_function = nn.MSELoss()
        scaler = torch.amp.GradScaler("cuda", enabled=training["amp"] and device.type == "cuda")
        model_dir.mkdir(parents=True, exist_ok=True)
        history = []
        best_loss, stale = math.inf, 0
        for epoch in range(1, training["epochs"] + 1):
            model.train(); total_loss = 0.0; count = 0
            for x, y in loaders["train"]:
                x = x.to(device, non_blocking=True); y = y.to(device, non_blocking=True)
                optimizer.zero_grad(set_to_none=True)
                with torch.autocast(device_type=device.type, dtype=torch.float16,
                                    enabled=training["amp"] and device.type == "cuda"):
                    prediction = model(x); loss = loss_function(prediction, y)
                scaler.scale(loss).backward(); scaler.step(optimizer); scaler.update()
                total_loss += float(loss.detach()) * len(x); count += len(x)
            val_y, val_p = _evaluate(model, loaders["validation"], device)
            val_loss = float(np.mean((val_y - val_p) ** 2))
            record = {"epoch": epoch, "train_mse": total_loss / count,
                      "validation_mse": val_loss, "seconds": int(time.time())}
            history.append(record)
            print(json.dumps({"species": name, **record}), flush=True)
            if val_loss < best_loss:
                best_loss, stale = val_loss, 0
                torch.save({"state_dict": model.state_dict(), "species": name,
                            "head_hidden": training["head_hidden"], "epoch": epoch,
                            "validation_mse": val_loss, "config": plan_model(config)}, best_path)
            else:
                stale += 1
                if stale >= training["patience"]:
                    break
        checkpoint = torch.load(best_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint["state_dict"])
        metrics = {}
        for split in ("validation", "test"):
            observed, predicted = _evaluate(model, loaders[split], device)
            metrics[split] = _metrics(observed, predicted)
            np.savez_compressed(model_dir / f"{split}_predictions.npz",
                                observed=observed, predicted=predicted)
        with (model_dir / "history.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=history[0], delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(history)
        metadata = {"species": name, "device": str(device), "torch": torch.__version__,
                    "best_epoch": checkpoint["epoch"], "best_validation_mse": best_loss,
                    "metrics": metrics, "training": training, "seed": settings["seed"],
                    "checkpoint_sha256": _sha256(best_path)}
        (model_dir / "metrics.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
        results[name] = metadata
    return results


def _iter_inference_windows(chrom: str, length: int, window_bp: int,
                            step_bp: int) -> Iterator[tuple[str, int, int]]:
    for start in range(0, length - window_bp + 1, step_bp):
        yield chrom, start, start + window_bp


def _iter_region_windows(chrom: str, length: int, window_bp: int, step_bp: int,
                         center_bp: int, region_start: int,
                         region_end: int) -> Iterator[tuple[str, int, int]]:
    flank = (window_bp - center_bp) // 2
    first = max(0, ((region_start - flank - center_bp) // step_bp) * step_bp)
    last = min(length - window_bp, math.ceil((region_end - flank) / step_bp) * step_bp)
    for start in range(first, last + 1, step_bp):
        output_start, output_end = start + flank, start + flank + center_bp
        if output_end > region_start and output_start < region_end:
            yield chrom, start, start + window_bp


def parse_region(value: str) -> tuple[str, int, int]:
    try:
        chrom, coordinates = value.rsplit(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        start, end = int(start_text.replace(",", "")), int(end_text.replace(",", ""))
    except (ValueError, AttributeError) as exc:
        raise ValueError("--region must have the form CHROM:START-END") from exc
    if not chrom or start < 0 or end <= start:
        raise ValueError("--region must have non-negative START and END > START")
    return chrom, start, end


def predict_species(config: dict, requested: Sequence[str] | None = None,
                    device_name: str | None = None, chromosomes: Sequence[str] | None = None,
                    force: bool = False, region: tuple[str, int, int] | None = None) -> dict:
    np, torch, _, _ = _lazy_model_imports()
    settings, training = config["settings"], config["training"]
    results = {}
    for name, item in select_species(config, requested):
        root = Path(config["output"]) / item["slug"]
        checkpoint_path = root / "model" / "best.pt"
        if not checkpoint_path.is_file():
            raise ModelWorkflowError(f"Missing checkpoint: {checkpoint_path}")
        prediction_dir = root / "predictions"
        summary_path = prediction_dir / "summary.json"
        if summary_path.exists() and not force:
            raise ModelWorkflowError(f"Refusing to overwrite {summary_path}; pass --force")
        prediction_dir.mkdir(parents=True, exist_ok=True)
        device = torch.device(device_name or ("cuda" if torch.cuda.is_available() else "cpu"))
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model = build_model(checkpoint["head_hidden"]).to(device)
        model.load_state_dict(checkpoint["state_dict"]); model.eval()
        fasta = IndexedFasta(item["fasta"])
        allowed = [chrom for split in ("train", "validation", "test")
                   for chrom in item["chromosomes"][split]]
        if region and chromosomes:
            raise ValueError("Use either --region or --chromosome, not both")
        selected_chroms = [region[0]] if region else (list(chromosomes) if chromosomes else allowed)
        unknown = set(selected_chroms) - set(allowed)
        if unknown:
            raise ModelWorkflowError(f"Chromosomes are outside configured nuclear set: {sorted(unknown)}")
        window_bp = settings["window_bp"]; step = settings["inference_step_bp"]
        center = settings["inference_center_bp"]; flank = (window_bp - center) // 2
        output_paths = {key: prediction_dir / f"{item['slug']}.{key}.bedGraph"
                        for key in ("plus", "minus", "minus_signed")}
        handles = {key: path.open("w", encoding="utf-8") for key, path in output_paths.items()}
        total_windows = 0
        skipped_ambiguous = 0
        try:
            batch_sequences = []; batch_coordinates = []
            def flush():
                nonlocal total_windows, batch_sequences, batch_coordinates
                if not batch_sequences:
                    return
                tensor = torch.from_numpy(np.stack(batch_sequences)).to(device)
                with torch.no_grad(), torch.autocast(device_type=device.type, dtype=torch.float16,
                                                     enabled=training["amp"] and device.type == "cuda"):
                    values = model(tensor).float().cpu().numpy()
                for (chrom, start, _), (plus, minus) in zip(batch_coordinates, values):
                    out_start, out_end = start + flank, start + flank + center
                    handles["plus"].write(f"{chrom}\t{out_start}\t{out_end}\t{plus:.7g}\n")
                    handles["minus"].write(f"{chrom}\t{out_start}\t{out_end}\t{minus:.7g}\n")
                    handles["minus_signed"].write(f"{chrom}\t{out_start}\t{out_end}\t{-minus:.7g}\n")
                total_windows += len(batch_sequences); batch_sequences = []; batch_coordinates = []
            for chrom in selected_chroms:
                length = fasta.index[chrom][0]
                coordinates = (_iter_region_windows(chrom, length, window_bp, step, center,
                                                     region[1], region[2]) if region else
                               _iter_inference_windows(chrom, length, window_bp, step))
                for coordinate in coordinates:
                    sequence = fasta.fetch(*coordinate)
                    if _ambiguous_fraction(sequence) > settings["max_ambiguous_fraction"]:
                        skipped_ambiguous += 1
                        continue
                    batch_sequences.append(one_hot(sequence)); batch_coordinates.append(coordinate)
                    if len(batch_sequences) >= max(training["batch_size"], 32):
                        flush()
                flush()
        finally:
            fasta.close()
            for handle in handles.values():
                handle.close()
        summary = {"species": name, "checkpoint": str(checkpoint_path),
                   "checkpoint_sha256": _sha256(checkpoint_path), "window_bp": window_bp,
                   "step_bp": step, "center_bp": center, "chromosomes": selected_chroms,
                   "region": ({"chrom": region[0], "start0": region[1], "end0": region[2]}
                              if region else None),
                   "windows": total_windows, "skipped_ambiguous_windows": skipped_ambiguous,
                   "outputs": {key: {"path": str(path), "sha256": _sha256(path)}
                               for key, path in output_paths.items()},
                   "value_scale": "raw sigmoid output; no threshold or min-max normalization"}
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        results[name] = summary
    return results


def _read_track(path: Path, chrom: str, start: int, end: int) -> tuple[list[float], list[float]]:
    x, y = [], []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != chrom:
                continue
            left, right = int(fields[1]), int(fields[2])
            if right <= start or left >= end:
                continue
            x.append((left + right) / 2); y.append(float(fields[3]))
    return x, y


def _read_observed(path: Path, chrom: str, start: int, end: int):
    x, plus, minus = [], [], []
    with _open_text(path) as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["chrom"] != chrom:
                continue
            left, right = int(row["start0"]), int(row["end0"])
            if right <= start or left >= end:
                continue
            x.append((left + right) / 2); plus.append(float(row["y_plus"])); minus.append(-float(row["y_minus"]))
    return x, plus, minus


def plot_region(config: dict, species_key: str, chrom: str, start: int, end: int,
                output: str | Path) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ModelWorkflowError("Matplotlib is required for plot") from exc
    selected = select_species(config, [species_key])
    name, item = selected[0]
    root = Path(config["output"]) / item["slug"]
    observed = _read_observed(Path(item["positive_y"]), chrom, start, end)
    plus = _read_track(root / "predictions" / f"{item['slug']}.plus.bedGraph", chrom, start, end)
    minus = _read_track(root / "predictions" / f"{item['slug']}.minus_signed.bedGraph", chrom, start, end)
    figure, axes = plt.subplots(2, 1, figsize=(12, 4.5), sharex=True, constrained_layout=True)
    axes[0].fill_between(observed[0], observed[1], 0, step="mid", color="#087f5b", alpha=0.85,
                         label="plus strand")
    axes[0].fill_between(observed[0], observed[2], 0, step="mid", color="#7048e8", alpha=0.85,
                         label="minus strand")
    axes[0].set_ylabel("Observed y"); axes[0].set_title(f"{name}  {chrom}:{start:,}-{end:,}")
    axes[1].fill_between(plus[0], plus[1], 0, step="mid", color="#20c997", alpha=0.8,
                         label="plus prediction")
    axes[1].fill_between(minus[0], minus[1], 0, step="mid", color="#9775fa", alpha=0.8,
                         label="minus prediction")
    axes[1].set_ylabel("Predicted y"); axes[1].set_xlabel(f"Genomic coordinate on {chrom} (bp)")
    for axis in axes:
        axis.axhline(0, color="black", linewidth=0.7); axis.set_xlim(start, end)
        limit = max(abs(value) for values in axis.collections for value in values.get_datalim(axis.transData).intervaly)
        axis.set_ylim(-max(1.0, limit), max(1.0, limit)); axis.legend(frameon=False, ncol=2, loc="upper right")
        axis.spines[["top", "right"]].set_visible(False)
    output_path = Path(output); output_path.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in (".png", ".pdf"):
        path = output_path.with_suffix(suffix); figure.savefig(path, dpi=300 if suffix == ".png" else None)
        paths.append(str(path))
    plt.close(figure)
    return paths


def _print_json(value) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train strand-resolved SeiPlant-style GRO-seq models")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("plan", "prepare", "train", "predict"):
        child = sub.add_parser(command)
        child.add_argument("config")
        child.add_argument("--species", action="append", help="Species name or slug; repeatable")
        if command in {"prepare", "train", "predict"}:
            child.add_argument("--force", action="store_true")
        if command in {"train", "predict"}:
            child.add_argument("--device", help="PyTorch device, e.g. cuda:0 or cpu")
        if command == "predict":
            child.add_argument("--chromosome", action="append")
            child.add_argument("--region", help="Optional bounded inference region, CHROM:START-END")
    plot = sub.add_parser("plot")
    plot.add_argument("config"); plot.add_argument("--species", required=True)
    plot.add_argument("--chrom", required=True); plot.add_argument("--start", type=int, required=True)
    plot.add_argument("--end", type=int, required=True); plot.add_argument("--output", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = load_model_config(args.config)
        if args.command == "plan":
            _print_json(plan_model(config))
        elif args.command == "prepare":
            _print_json(prepare_manifests(config, args.species, args.force))
        elif args.command == "train":
            _print_json(train_species(config, args.species, args.device, args.force))
        elif args.command == "predict":
            region = parse_region(args.region) if args.region else None
            _print_json(predict_species(config, args.species, args.device, args.chromosome,
                                        args.force, region))
        elif args.command == "plot":
            if args.end <= args.start:
                raise ValueError("--end must be greater than --start")
            _print_json({"outputs": plot_region(config, args.species, args.chrom,
                                                  args.start, args.end, args.output)})
        return 0
    except (ValueError, OSError, ModelWorkflowError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
