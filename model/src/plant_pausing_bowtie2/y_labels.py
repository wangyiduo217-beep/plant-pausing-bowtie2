"""Reproducible strand-resolved GRO-seq interval and y-label generation.

The implementation mirrors the production method documented in
docs/groseq-y-label-generation.md while keeping paths and libraries in JSON.
Execution requires Linux, HOMER, SAMtools and BEDTools. Planning and unit tests
do not invoke those tools or read large BAM files.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import concurrent.futures
import copy
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid


VERSION = "groseq_y_article_analog_v2.0"
HOMER_PARAMETERS = [
    "-style", "groseq", "-tssFold", "4", "-bodyFold", "3",
    "-minBodySize", "500", "-maxBodySize", "100000",
    "-pseudoCount", "1", "-o", "auto",
]
DEFAULT_TOOLS = {
    "samtools": "samtools",
    "bedtools": "bedtools",
    "make_tag_directory": "makeTagDirectory",
    "find_peaks": "findPeaks",
    "pos2bed": "pos2bed.pl",
}
SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")


class LabelWorkflowError(RuntimeError):
    """Raised when inputs, external commands or output checks fail."""


def _utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _resolve(value: object, base: Path, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty path string")
    path = Path(value).expanduser()
    return str((base / path).resolve() if not path.is_absolute() else path.resolve())


def _positive_int(value: object, name: str, minimum: int = 1) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def load_y_config(path: str | os.PathLike[str]) -> dict:
    """Load and validate a y-label config; relative paths use its directory."""
    config_path = Path(path).expanduser().resolve()
    raw = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Configuration root must be an object")
    base = config_path.parent
    output = _resolve(raw.get("output", "../results/groseq_y"), base, "output")

    settings = copy.deepcopy(raw.get("settings", {}))
    if not isinstance(settings, dict):
        raise ValueError("settings must be an object")
    allowed_settings = {
        "window_bp", "step_bp", "minimum_support", "threads", "jobs",
        "keep_intermediates", "tools",
    }
    unknown = set(settings) - allowed_settings
    if unknown:
        raise ValueError("Unsupported settings: " + ", ".join(sorted(unknown)))
    normalized_settings = {
        "window_bp": _positive_int(settings.get("window_bp", 1024), "window_bp"),
        "step_bp": _positive_int(settings.get("step_bp", 512), "step_bp"),
        "minimum_support": _positive_int(settings.get("minimum_support", 2), "minimum_support", 2),
        "threads": _positive_int(settings.get("threads", 8), "threads"),
        "jobs": _positive_int(settings.get("jobs", 3), "jobs"),
        "keep_intermediates": settings.get("keep_intermediates", False),
    }
    if type(normalized_settings["keep_intermediates"]) is not bool:
        raise ValueError("keep_intermediates must be boolean")
    tools = dict(DEFAULT_TOOLS)
    supplied_tools = settings.get("tools", {})
    if not isinstance(supplied_tools, dict) or set(supplied_tools) - set(tools):
        raise ValueError("tools contains an unsupported key")
    for key, value in supplied_tools.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"tools.{key} must be a nonempty executable name or path")
        tools[key] = value
    normalized_settings["tools"] = tools

    species_raw = raw.get("species")
    if not isinstance(species_raw, dict) or not species_raw:
        raise ValueError("species must be a nonempty object")
    species: dict[str, dict] = {}
    slugs: set[str] = set()
    for name, entry in species_raw.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(entry, dict):
            raise ValueError("Every species entry must have a nonempty name and object value")
        slug = entry.get("slug")
        if not isinstance(slug, str) or not SAFE_ID.fullmatch(slug) or slug in slugs:
            raise ValueError(f"Invalid or duplicate species slug for {name}")
        slugs.add(slug)
        status = entry.get("dataset_status", "article_analog_primary")
        if not isinstance(status, str) or not status.strip() or any(ord(c) < 32 for c in status):
            raise ValueError(f"Invalid dataset_status for {name}")
        mask = entry.get("mask")
        species[name] = {
            "slug": slug,
            "fai": _resolve(entry.get("fai"), base, f"species.{name}.fai"),
            "mask": _resolve(mask, base, f"species.{name}.mask") if mask else None,
            "dataset_status": status,
        }

    libraries_raw = raw.get("libraries")
    if not isinstance(libraries_raw, list) or not libraries_raw:
        raise ValueError("libraries must be a nonempty list")
    libraries = []
    seen: set[str] = set()
    for entry in libraries_raw:
        if not isinstance(entry, dict):
            raise ValueError("Each library must be an object")
        srx = entry.get("srx")
        if not isinstance(srx, str) or not SAFE_ID.fullmatch(srx) or srx in seen:
            raise ValueError("Every library must have a distinct safe SRX identifier")
        seen.add(srx)
        species_name = entry.get("species")
        if species_name not in species:
            raise ValueError(f"Unknown species for {srx}: {species_name}")
        layout = entry.get("layout")
        if layout not in ("SINGLE", "PAIRED"):
            raise ValueError(f"layout for {srx} must be SINGLE or PAIRED")
        bams = entry.get("bams")
        runs = entry.get("runs")
        if not isinstance(bams, list) or not bams or any(not isinstance(x, str) for x in bams):
            raise ValueError(f"bams for {srx} must be a nonempty path list")
        if not isinstance(runs, list) or len(runs) != len(bams) or any(
                not isinstance(x, str) or not SAFE_ID.fullmatch(x) for x in runs):
            raise ValueError(f"runs for {srx} must contain one safe identifier per BAM")
        reverse = entry.get("reverse_strand", False)
        if type(reverse) is not bool:
            raise ValueError(f"reverse_strand for {srx} must be boolean")
        libraries.append({
            "srx": srx,
            "species": species_name,
            "project": entry.get("project"),
            "group": entry.get("group"),
            "layout": layout,
            "runs": runs,
            "bams": [_resolve(x, base, f"bams for {srx}") for x in bams],
            "reverse_strand": reverse,
        })
    counts = collections.Counter(x["species"] for x in libraries)
    for name in species:
        if counts[name] < normalized_settings["minimum_support"]:
            raise ValueError(
                f"{name} has {counts[name]} libraries, below minimum_support "
                f"{normalized_settings['minimum_support']}"
            )
    return {
        "config_path": str(config_path),
        "output": output,
        "settings": normalized_settings,
        "species": species,
        "libraries": libraries,
    }


def _library_commands(config: dict, library: dict) -> dict:
    settings, tools = config["settings"], config["settings"]["tools"]
    root = Path(config["output"])
    mask = config["species"][library["species"]]["mask"]
    analysis = []
    for run, bam in zip(library["runs"], library["bams"]):
        target = root / "work" / library["srx"] / f"{run}.analysis_input.bam"
        view = [tools["samtools"], "view", "-@", str(settings["threads"]), "-h", "-b"]
        if library["layout"] == "PAIRED":
            view += ["-f", "64"]
        view.append(bam)
        if mask:
            command = {
                "pipeline": [
                    view,
                    [tools["bedtools"], "intersect", "-v", "-abam", "stdin", "-b", mask],
                ],
                "stdout": str(target),
            }
        else:
            command = {"argv": view + ["-o", str(target)]}
        analysis.append({"run": run, "bam": bam, "output": str(target), "command": command})
    tag_dir = root / "work" / library["srx"] / "tag_directory"
    make_tags = [tools["make_tag_directory"], str(tag_dir)] + [x["output"] for x in analysis]
    find_peaks = [tools["find_peaks"], str(tag_dir)] + HOMER_PARAMETERS
    if library["reverse_strand"]:
        find_peaks.append("-rev")
    return {
        "srx": library["srx"],
        "analysis_inputs": analysis,
        "make_tag_directory": make_tags,
        "find_peaks": find_peaks,
        "pos2bed": [tools["pos2bed"], str(tag_dir / "transcripts.txt")],
        "output": str(root / "libraries" / library["srx"] / f"{library['srx']}.transcripts.bed"),
    }


def plan_y_labels(config: dict) -> dict:
    """Return the external commands and scientific settings without execution."""
    return {
        "version": VERSION,
        "output": config["output"],
        "homer_parameters": HOMER_PARAMETERS,
        "window_bp": config["settings"]["window_bp"],
        "step_bp": config["settings"]["step_bp"],
        "minimum_support": config["settings"]["minimum_support"],
        "libraries": [_library_commands(config, x) for x in config["libraries"]],
        "build_y": [
            {
                "species": name,
                "slug": entry["slug"],
                "independent_libraries": sum(x["species"] == name for x in config["libraries"]),
                "strand_operations": "split -> bedtools sort -> bedtools multiinter; plus and minus separately",
                "fai": entry["fai"],
                "mask": entry["mask"],
                "dataset_status": entry["dataset_status"],
            }
            for name, entry in config["species"].items()
        ],
    }


def _resolve_tools(config: dict, required: set[str]) -> dict[str, str]:
    resolved = {}
    for key in required:
        value = config["settings"]["tools"][key]
        executable = shutil.which(value)
        if executable is None and Path(value).is_file():
            executable = str(Path(value).resolve())
        if executable is None:
            raise LabelWorkflowError(f"Required executable is missing: {value}")
        resolved[key] = executable
    return resolved


def _identity(path: str | Path) -> dict:
    item = Path(path).resolve(strict=True)
    if not item.is_file() or item.stat().st_size == 0:
        raise LabelWorkflowError(f"Expected a nonempty file: {item}")
    stat = item.stat()
    return {"path": str(item), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _fingerprint(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def _run(argv: list[str], log, stdout_path: Path | None = None) -> None:
    log.write(_utc() + " COMMAND " + shlex.join(list(map(str, argv))) + "\n")
    log.flush()
    output = stdout_path.open("wb") if stdout_path else None
    try:
        subprocess.run(list(map(str, argv)), check=True, stdout=output or log,
                       stderr=log, stdin=subprocess.DEVNULL)
    finally:
        if output:
            output.close()


def _pipeline(first: list[str], second: list[str], output: Path, log) -> None:
    log.write(_utc() + " PIPE " + shlex.join(first) + " | " + shlex.join(second) + "\n")
    log.flush()
    left = subprocess.Popen(first, stdout=subprocess.PIPE, stderr=log, stdin=subprocess.DEVNULL)
    try:
        with output.open("wb") as stream:
            right = subprocess.run(second, stdin=left.stdout, stdout=stream, stderr=log)
        left.stdout.close()
        left_code = left.wait()
    except BaseException:
        left.kill()
        left.wait()
        raise
    if left_code or right.returncode:
        raise subprocess.CalledProcessError(left_code or right.returncode, first if left_code else second)


def _mask_is_nonempty(path: str | None) -> bool:
    return bool(path and Path(path).is_file() and Path(path).stat().st_size > 0)


def _call_one(config: dict, library: dict, tools: dict[str, str]) -> dict:
    output = Path(config["output"])
    species = config["species"][library["species"]]
    final = output / "libraries" / library["srx"] / f"{library['srx']}.transcripts.bed"
    metadata_path = output / "metadata" / f"{library['srx']}.json"
    log_path = output / "logs" / f"{library['srx']}.log"
    bam_identities = [_identity(x) for x in library["bams"]]
    for bam in library["bams"]:
        index = Path(bam + ".csi")
        if not index.is_file() or index.stat().st_size == 0:
            raise LabelWorkflowError(f"CSI index is missing for {bam}")
    mask_identity = _identity(species["mask"]) if _mask_is_nonempty(species["mask"]) else None
    inputs = {
        "library": library,
        "bams": bam_identities,
        "mask": mask_identity,
        "homer_parameters": HOMER_PARAMETERS,
        "tools": {name: _identity(path) for name, path in sorted(tools.items())},
    }
    fingerprint = _fingerprint(inputs)
    if metadata_path.is_file():
        old = json.loads(metadata_path.read_text(encoding="utf-8"))
        if old.get("state") == "complete" and old.get("fingerprint") == fingerprint and final.is_file():
            return old
        raise LabelWorkflowError(f"Existing output for {library['srx']} does not match this configuration")
    if final.exists():
        raise LabelWorkflowError(f"Untracked output already exists: {final}")

    token = uuid.uuid4().hex
    work = output / "work" / f"{library['srx']}.building-{token}"
    tag_dir = work / "tag_directory"
    work.mkdir(parents=True, exist_ok=False)
    final.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    analysis_bams = []
    try:
        with log_path.open("a", encoding="utf-8") as log:
            for run, source in zip(library["runs"], library["bams"]):
                target = work / f"{run}.analysis_input.bam"
                view = [tools["samtools"], "view", "-@", str(config["settings"]["threads"]), "-h", "-b"]
                if library["layout"] == "PAIRED":
                    view += ["-f", "64"]
                view.append(source)
                if _mask_is_nonempty(species["mask"]):
                    intersect = [tools["bedtools"], "intersect", "-v", "-abam", "stdin", "-b", species["mask"]]
                    _pipeline(view, intersect, target, log)
                else:
                    _run(view + ["-o", str(target)], log)
                _run([tools["samtools"], "quickcheck", "-v", str(target)], log)
                analysis_bams.append(target)

            _run([tools["make_tag_directory"], str(tag_dir), *map(str, analysis_bams)], log)
            find = [tools["find_peaks"], str(tag_dir), *HOMER_PARAMETERS]
            if library["reverse_strand"]:
                find.append("-rev")
            _run(find, log)
            transcripts = tag_dir / "transcripts.txt"
            if not transcripts.is_file() or transcripts.stat().st_size == 0:
                raise LabelWorkflowError(f"HOMER did not create transcripts.txt for {library['srx']}")
            converted = work / "converted.bed"
            _run([tools["pos2bed"], str(transcripts)], log, converted)
            filtered = work / "filtered.bed"
            if _mask_is_nonempty(species["mask"]):
                _run([tools["bedtools"], "intersect", "-v", "-a", str(converted),
                      "-b", species["mask"]], log, filtered)
            else:
                shutil.copyfile(converted, filtered)

        strands = {"+": 0, "-": 0}
        rows = 0
        with filtered.open(encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[5] not in strands:
                    raise LabelWorkflowError(f"Invalid BED row for {library['srx']}: {line[:160]}")
                int(fields[1]); int(fields[2])
                if int(fields[2]) <= int(fields[1]):
                    raise LabelWorkflowError(f"Invalid BED coordinates for {library['srx']}")
                strands[fields[5]] += 1
                rows += 1
        if rows == 0 or min(strands.values()) == 0:
            raise LabelWorkflowError(f"Empty or single-strand interval output for {library['srx']}")
        installing = final.with_name(final.name + ".installing")
        shutil.copyfile(filtered, installing)
        installing.replace(final)
        result = {
            "version": VERSION,
            "state": "complete",
            "completed_utc": _utc(),
            "fingerprint": fingerprint,
            "srx": library["srx"],
            "species": library["species"],
            "project": library.get("project"),
            "group": library.get("group"),
            "runs": library["runs"],
            "layout": library["layout"],
            "read_selection": "R1 only (SAM flag 64)" if library["layout"] == "PAIRED" else "all primary single-end records",
            "strand_interpretation": "reverse alignment strand (-rev)" if library["reverse_strand"] else "alignment strand equals RNA strand",
            "homer_parameters": HOMER_PARAMETERS,
            "inputs": inputs,
            "intervals": rows,
            "by_strand": strands,
            "output": str(final),
            "log": str(log_path),
        }
        _write_json(metadata_path, result)
        if not config["settings"]["keep_intermediates"]:
            shutil.rmtree(work)
        return result
    except BaseException:
        # Preserve the isolated work directory and log for diagnosis.
        raise


def call_intervals(config: dict) -> list[dict]:
    if not sys.platform.startswith("linux"):
        raise LabelWorkflowError("Interval execution requires Linux; use plan on other systems")
    tools = _resolve_tools(config, set(DEFAULT_TOOLS))
    Path(config["output"]).mkdir(parents=True, exist_ok=True)
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=config["settings"]["jobs"]) as pool:
        futures = {pool.submit(_call_one, config, library, tools): library["srx"]
                   for library in config["libraries"]}
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda x: x["srx"])
    _write_json(Path(config["output"]) / "intervals_summary.json", {
        "version": VERSION,
        "state": "complete",
        "completed_utc": _utc(),
        "libraries": results,
    })
    return results


def confidence_for_support(support: int, observed_supports: list[int],
                           total_libraries: int | None = None) -> float:
    values = sorted(set(observed_supports))
    if not values or support not in values:
        raise ValueError("support must occur in observed_supports")
    if len(values) == 1:
        return 1.0
    if total_libraries is None:
        total_libraries = max(values)
    _positive_int(total_libraries, "total_libraries")
    # Keep the production operation order. Although total_libraries cancels
    # algebraically, simplifying this expression changes binary floating-point
    # values at a few rounding boundaries in large real data sets.
    fraction = support / total_libraries
    minimum = values[0] / total_libraries
    maximum = values[-1] / total_libraries
    return 0.1 + 0.9 * (fraction - minimum) / (maximum - minimum)


def rounded_y(raw: float) -> float:
    if raw < 0:
        raise ValueError("raw y must be nonnegative")
    value = min(math.floor(raw * 10 + 0.5) / 10, 1.0)
    return 0.0 if value < 0.1 else value


def _read_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise LabelWorkflowError(f"Invalid FAI row: {line[:160]}")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise LabelWorkflowError(f"Empty FAI: {path}")
    return lengths


def _read_merged_masks(path: str | None) -> dict[str, list[tuple[int, int]]]:
    masks: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    if not _mask_is_nonempty(path):
        return masks
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise LabelWorkflowError(f"Invalid mask row: {line[:160]}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise LabelWorkflowError(f"Invalid mask coordinates: {line[:160]}")
            masks[fields[0]].append((start, end))
    for chrom, rows in masks.items():
        merged = []
        for start, end in sorted(rows):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        masks[chrom] = merged
    return masks


def _overlaps_mask(chrom: str, start: int, end: int,
                   masks: dict[str, list[tuple[int, int]]], starts: dict[str, list[int]]) -> bool:
    rows = masks.get(chrom, [])
    index = bisect.bisect_left(starts.get(chrom, []), end)
    return index > 0 and rows[index - 1][1] > start


def score_windows(intervals: list[dict], lengths: dict[str, int], window: int,
                  step: int) -> dict[tuple[str, int, int], list[float | int]]:
    """Pure sliding-window implementation used by production and unit tests."""
    scores: dict[tuple[str, int, int], list[float | int]] = collections.defaultdict(
        lambda: [0.0, 0.0, 0, 0]
    )
    for item in intervals:
        chrom, left, right = item["chrom"], item["start"], item["end"]
        if chrom not in lengths:
            continue
        first = max(0, (left - window) // step + 1)
        last = (right - 1) // step
        axis = 0 if item["strand"] == "+" else 1
        for index in range(first, last + 1):
            start, end = index * step, index * step + window
            if end > lengths[chrom]:
                continue
            overlap = max(0, min(end, right) - max(start, left))
            if overlap == 0:
                continue
            key = (chrom, start, end)
            scores[key][axis] += overlap / window * item["confidence"]
            scores[key][2 + axis] = max(int(scores[key][2 + axis]), item["support"])
    return scores


def _species_intervals(config: dict, species_name: str, stage: Path,
                       bedtools: str, log) -> list[dict]:
    libraries = [x for x in config["libraries"] if x["species"] == species_name]
    output = Path(config["output"])
    strand_files: dict[str, list[tuple[str, Path]]] = {"+": [], "-": []}
    for library in libraries:
        source = output / "libraries" / library["srx"] / f"{library['srx']}.transcripts.bed"
        if not source.is_file() or source.stat().st_size == 0:
            raise LabelWorkflowError(f"Missing interval BED for {library['srx']}: {source}")
        handles = {strand: (stage / f"{library['srx']}.{strand}.unsorted.bed").open("w", encoding="utf-8")
                   for strand in ("+", "-")}
        try:
            with source.open(encoding="utf-8") as stream:
                for line in stream:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 6 or fields[5] not in handles:
                        raise LabelWorkflowError(f"Invalid interval row in {source}: {line[:160]}")
                    handles[fields[5]].write(line)
        finally:
            for handle in handles.values():
                handle.close()
        for strand in ("+", "-"):
            unsorted = stage / f"{library['srx']}.{strand}.unsorted.bed"
            sorted_path = stage / f"{library['srx']}.{strand}.bed"
            _run([bedtools, "sort", "-i", str(unsorted)], log, sorted_path)
            strand_files[strand].append((library["srx"], sorted_path))

    intervals = []
    minimum = config["settings"]["minimum_support"]
    for strand in ("+", "-"):
        raw = stage / f"multiinter.{strand}.tsv"
        command = [bedtools, "multiinter", "-i", *[str(x[1]) for x in strand_files[strand]],
                   "-names", *[x[0] for x in strand_files[strand]]]
        _run(command, log, raw)
        retained = []
        with raw.open(encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise LabelWorkflowError(f"Invalid multiinter row: {line[:160]}")
                support = int(fields[3])
                if support >= minimum:
                    retained.append((fields[0], int(fields[1]), int(fields[2]), support, fields[4]))
        observed = [x[3] for x in retained]
        for chrom, start, end, support, names in retained:
            intervals.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "support": support,
                "supporting_srx": names,
                "confidence": confidence_for_support(support, observed, len(libraries)),
            })
    intervals.sort(key=lambda x: (x["chrom"], x["start"], x["end"], x["strand"]))
    if not intervals:
        raise LabelWorkflowError(f"No consensus intervals passed minimum support for {species_name}")
    return intervals


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        block = stream.read(1024 * 1024)
        while block:
            digest.update(block)
            block = stream.read(1024 * 1024)
    return digest.hexdigest()


def build_y(config: dict) -> dict:
    if not sys.platform.startswith("linux"):
        raise LabelWorkflowError("y-label execution requires Linux; use plan on other systems")
    tools = _resolve_tools(config, {"bedtools"})
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=True)
    library_outputs = []
    for library in config["libraries"]:
        path = output / "libraries" / library["srx"] / f"{library['srx']}.transcripts.bed"
        library_outputs.append(_identity(path))
    fais = {name: _identity(entry["fai"]) for name, entry in config["species"].items()}
    masks = {name: _identity(entry["mask"]) if _mask_is_nonempty(entry["mask"]) else None
             for name, entry in config["species"].items()}
    build_identity = {
        "version": VERSION,
        "settings": {key: config["settings"][key] for key in ("window_bp", "step_bp", "minimum_support")},
        "species": config["species"],
        "libraries": config["libraries"],
        "interval_files": library_outputs,
        "fai_files": fais,
        "mask_files": masks,
        "bedtools": _identity(tools["bedtools"]),
    }
    fingerprint = _fingerprint(build_identity)
    summary_path = output / "groseq_y_summary.json"
    if summary_path.is_file():
        old = json.loads(summary_path.read_text(encoding="utf-8"))
        if old.get("state") == "complete" and old.get("fingerprint") == fingerprint:
            return old
        raise LabelWorkflowError("Existing y release does not match this configuration; use a new output directory")
    labels_target = output / "labels"
    if labels_target.exists():
        raise LabelWorkflowError(f"Untracked labels directory already exists: {labels_target}")

    stage = output / f".labels-building-{uuid.uuid4().hex}"
    stage_labels = stage / "labels"
    stage_labels.mkdir(parents=True)
    log_path = output / "logs" / "build_y.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "version": VERSION,
        "state": "building",
        "fingerprint": fingerprint,
        "created_utc": _utc(),
        "window_bp": config["settings"]["window_bp"],
        "step_bp": config["settings"]["step_bp"],
        "minimum_independent_SRX_support": config["settings"]["minimum_support"],
        "interval_caller": "HOMER findPeaks -style groseq",
        "homer_parameters": HOMER_PARAMETERS,
        "species": {},
        "interpretation": "strand-specific reproducible nascent-transcription interval labels; not pausing index or single-base pause truth",
    }
    try:
        with log_path.open("a", encoding="utf-8") as log:
            for species_name, species in config["species"].items():
                species_stage = stage / species["slug"]
                species_stage.mkdir()
                intervals = _species_intervals(config, species_name, species_stage, tools["bedtools"], log)
                consensus = stage_labels / f"{species['slug']}.consensus_intervals.bed"
                with consensus.open("w", encoding="utf-8") as stream:
                    for index, item in enumerate(intervals, 1):
                        stream.write(
                            f"{item['chrom']}\t{item['start']}\t{item['end']}\t"
                            f"{species['slug']}_consensus_{index}\t{item['confidence']:.6f}\t"
                            f"{item['strand']}\t{item['support']}\t{item['supporting_srx']}\n"
                        )
                lengths = _read_lengths(Path(species["fai"]))
                mask_rows = _read_merged_masks(species["mask"])
                mask_starts = {chrom: [x[0] for x in rows] for chrom, rows in mask_rows.items()}
                scores = score_windows(
                    intervals, lengths, config["settings"]["window_bp"], config["settings"]["step_bp"]
                )
                y_path = stage_labels / f"{species['slug']}.positive_windows.y.tsv.gz"
                kept = excluded = 0
                with gzip.open(y_path, "wt", encoding="utf-8") as stream:
                    stream.write(
                        "species\tchrom\tstart0\tend0\ty_plus\ty_minus\tbinary_plus\t"
                        "binary_minus\tmax_support_plus\tmax_support_minus\tdataset_status\n"
                    )
                    for (chrom, start, end), values in sorted(scores.items(), key=lambda x: (x[0][0], x[0][1])):
                        if _overlaps_mask(chrom, start, end, mask_rows, mask_starts):
                            excluded += 1
                            continue
                        plus, minus = rounded_y(float(values[0])), rounded_y(float(values[1]))
                        if plus == 0 and minus == 0:
                            continue
                        stream.write(
                            f"{species_name}\t{chrom}\t{start}\t{end}\t{plus:.1f}\t{minus:.1f}\t"
                            f"{int(float(values[0]) > 0)}\t{int(float(values[1]) > 0)}\t"
                            f"{int(values[2])}\t{int(values[3])}\t{species['dataset_status']}\n"
                        )
                        kept += 1
                selected = [x for x in config["libraries"] if x["species"] == species_name]
                summary["species"][species_name] = {
                    "libraries": [x["srx"] for x in selected],
                    "independent_libraries": len(selected),
                    "runs": sum(len(x["runs"]) for x in selected),
                    "consensus_atomic_intervals": len(intervals),
                    "positive_windows": kept,
                    "windows_excluded_by_mask": excluded,
                    "consensus_bed": str(labels_target / consensus.name),
                    "y_table": str(labels_target / y_path.name),
                    "dataset_status": species["dataset_status"],
                }
        labels_target.parent.mkdir(parents=True, exist_ok=True)
        stage_labels.replace(labels_target)
        summary["state"] = "complete"
        summary["completed_utc"] = _utc()
        _write_json(summary_path, summary)
        checksum_targets = sorted(labels_target.iterdir()) + [summary_path]
        with (output / "SHA256SUMS").open("w", encoding="utf-8") as stream:
            for path in checksum_targets:
                stream.write(f"{_sha256(path)}  {path.relative_to(output)}\n")
        shutil.rmtree(stage)
        return summary
    except BaseException:
        # Preserve the staging directory and log for inspection.
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Call stranded GRO-seq intervals and build 1024/512 y labels")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command, help_text in (
        ("plan", "Show commands and settings without reading BAM files"),
        ("call", "Call per-SRX stranded GRO-seq intervals"),
        ("build", "Build cross-library consensus intervals and y tables"),
        ("run", "Run interval calling followed by y-table generation"),
    ):
        sub = subparsers.add_parser(command, help=help_text)
        sub.add_argument("config", help="JSON y-label configuration")
    args = parser.parse_args(argv)
    try:
        config = load_y_config(args.config)
        if args.command == "plan":
            result = plan_y_labels(config)
        elif args.command == "call":
            result = {"libraries": call_intervals(config)}
        elif args.command == "build":
            result = build_y(config)
        else:
            call_intervals(config)
            result = build_y(config)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except (ValueError, OSError, LabelWorkflowError, subprocess.SubprocessError, KeyError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
