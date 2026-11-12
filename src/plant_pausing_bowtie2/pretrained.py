"""Verified pretrained downloads and explicit 1,024-bp FASTA prediction.

The original strand_model architecture and A/T/C/G encoding are reused unchanged.
Download and input validation require only the Python standard library.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

SPECIES = {
    "arabidopsis_thaliana": "Arabidopsis thaliana",
    "triticum_aestivum": "Triticum aestivum",
    "zea_mays": "Zea mays",
}
FIELDS = ("sequence_id", "prediction_plus", "prediction_minus")
WINDOW_BP = 1024
DOWNLOAD_CHUNK_BYTES = 4 * 1024 * 1024
DOWNLOAD_ATTEMPTS = 3


class PretrainedError(RuntimeError):
    """An input, checkpoint, download, or expected output failed validation."""


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("schema_version") != 1 or manifest.get("model_version") != "strand_seiplant_v1":
        raise PretrainedError("Unsupported pretrained manifest version or model identity")
    models = manifest.get("models")
    if not isinstance(models, dict) or set(models) != set(SPECIES):
        raise PretrainedError("Manifest must identify all three publication species")
    for slug, item in models.items():
        if not isinstance(item, dict):
            raise PretrainedError("Invalid model entry: " + slug)
        if not re.fullmatch(r"[0-9a-f]{64}", str(item.get("sha256", ""))):
            raise PretrainedError("Invalid checkpoint SHA-256: " + slug)
        if type(item.get("bytes")) is not int or item["bytes"] <= 0:
            raise PretrainedError("Invalid checkpoint byte count: " + slug)
        url = urllib.parse.urlsplit(str(item.get("url", "")))
        if url.scheme != "https" or not url.hostname or url.username or url.password:
            raise PretrainedError("Checkpoint downloads must use a public HTTPS URL")
        expected_srx = {"arabidopsis_thaliana": 14, "triticum_aestivum": 5, "zea_mays": 3}[slug]
        if item.get("independent_SRX") != expected_srx:
            raise PretrainedError("Wrong SRX cohort for publication checkpoint: " + slug)
        # Check a relative path here; checkpoint_path also checks symlink escapes.
        relative = Path(str(item.get("checkpoint_relative_path", "")))
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise PretrainedError("Checkpoint path must stay inside the repository")
    return manifest


def checkpoint_path(root: Path, item: dict) -> Path:
    root = root.resolve()
    target = (root / item["checkpoint_relative_path"]).resolve()
    try:
        target.relative_to(root)
    except ValueError as exc:
        raise PretrainedError("Checkpoint path escapes the repository") from exc
    if target == root:
        raise PretrainedError("Checkpoint path must name a file")
    return target


def verify_checkpoint(path: Path, item: dict) -> None:
    if not path.is_file():
        raise PretrainedError("Missing checkpoint: {}. Run scripts/download_pretrained.py first.".format(path))
    if path.stat().st_size != item["bytes"] or sha256_file(path) != item["sha256"]:
        raise PretrainedError("Checkpoint size or SHA-256 differs from the published model: " + str(path))


def _download_to_temporary(output, item: dict) -> str:
    """Bounded range requests prevent proxy truncation of long responses.

    Only complete, correctly positioned ranges are appended. A server that
    ignores the initial Range header may return a single full response. Once
    ranges have begun, a full response is rejected rather than concatenated.
    Every transfer still requires the final published byte count and SHA-256.
    """
    offset = 0
    url = item["url"]
    while offset < item["bytes"]:
        end = min(offset + DOWNLOAD_CHUNK_BYTES, item["bytes"]) - 1
        for attempt in range(DOWNLOAD_ATTEMPTS):
            output.seek(offset)
            output.truncate(offset)
            request = urllib.request.Request(url, headers={
                "User-Agent": "plant-pausing-bowtie2-pretrained",
                "Accept-Encoding": "identity",
                "Range": "bytes={}-{}".format(offset, end),
            })
            try:
                with urllib.request.urlopen(request, timeout=60) as response:
                    status = getattr(response, "status", 200)
                    if status == 206:
                        expected = "bytes {}-{}/{}".format(offset, end, item["bytes"])
                        if response.headers.get("Content-Range") != expected:
                            raise PretrainedError("Download returned an incorrect Content-Range")
                        block = response.read(end - offset + 2)
                        expected_size = end - offset + 1
                        if len(block) < expected_size:
                            raise http.client.IncompleteRead(b"", expected_size - len(block))
                        if len(block) != expected_size:
                            raise PretrainedError("Download range exceeds its declared byte count")
                        output.write(block)
                        offset = end + 1
                        # Reuse a verified HTTPS CDN redirect while it remains valid.
                        final_url = response.geturl()
                        if urllib.parse.urlsplit(final_url).scheme == "https":
                            url = final_url
                        break
                    if status != 200 or offset != 0:
                        raise PretrainedError("Server stopped honoring byte ranges or returned an unexpected status")
                    count = 0
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        count += len(block)
                        if count > item["bytes"]:
                            raise PretrainedError("Download exceeds the published checkpoint size")
                        output.write(block)
                    if count != item["bytes"]:
                        raise http.client.IncompleteRead(b"", item["bytes"] - count)
                    return "full_response"
            except urllib.error.HTTPError as exc:
                # Refresh an expired signed CDN URL; never publish partial bytes.
                if exc.code == 403 and url != item["url"]:
                    url = item["url"]
                elif exc.code != 429 and not 500 <= exc.code <= 599:
                    raise PretrainedError("Checkpoint request failed with HTTP " + str(exc.code)) from exc
                if attempt + 1 == DOWNLOAD_ATTEMPTS:
                    raise PretrainedError("Checkpoint request failed after three attempts") from exc
                time.sleep(0.5 * (attempt + 1))
            except (OSError, http.client.IncompleteRead) as exc:
                if attempt + 1 == DOWNLOAD_ATTEMPTS:
                    raise PretrainedError("Checkpoint transfer was interrupted after three attempts") from exc
                time.sleep(0.5 * (attempt + 1))
        else:
            raise PretrainedError("Checkpoint transfer exhausted its retry limit")
    return "verified_http_ranges"


def download_checkpoint(root: Path, item: dict, force: bool = False) -> dict:
    """Download with bounded retries; publish locally only after size/hash checks."""
    target = checkpoint_path(root, item)
    if target.exists():
        try:
            verify_checkpoint(target, item)
        except PretrainedError:
            if not force:
                raise PretrainedError("Existing checkpoint differs; inspect it or pass --force: " + str(target))
        else:
            return {"path": str(target), "sha256": item["sha256"], "status": "already_verified"}
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".partial", dir=str(target.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as output:
            transport = _download_to_temporary(output, item)
        if temporary.stat().st_size != item["bytes"] or sha256_file(temporary) != item["sha256"]:
            raise PretrainedError("Downloaded checkpoint failed size/SHA-256 verification")
        os.replace(str(temporary), str(target))
    finally:
        if temporary.exists():
            temporary.unlink()
    return {"path": str(target), "sha256": item["sha256"], "status": "downloaded_and_verified", "transport": transport}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read forward-reference windows; reject ambiguity, duplicate IDs and gaps."""
    records = []
    seen = set()
    identifier = None
    sequence_parts = []

    def finish():
        if identifier is None:
            return
        sequence = "".join(sequence_parts).upper()
        if len(sequence) != WINDOW_BP:
            raise PretrainedError("{}: expected 1024 bp, found {}".format(identifier, len(sequence)))
        if set(sequence) - set("ACGT"):
            raise PretrainedError(identifier + ": only A/C/G/T bases are accepted; N and gaps are rejected")
        records.append((identifier, sequence))

    with path.open(encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                tokens = line[1:].split()
                if not tokens:
                    raise PretrainedError("Empty FASTA identifier at line " + str(line_number))
                identifier = tokens[0]
                if identifier in seen:
                    raise PretrainedError("Duplicate FASTA identifier: " + identifier)
                seen.add(identifier)
                sequence_parts = []
            elif identifier is None:
                raise PretrainedError("Sequence encountered before a FASTA header")
            else:
                sequence_parts.append(line)
    finish()
    if not records:
        raise PretrainedError("FASTA has no sequences")
    return records


def _write_predictions(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def predict_fasta(*, root: Path, manifest: dict, species: str, fasta: Path,
                  output: Path, checkpoint: Path | None = None, device_name: str = "cpu",
                  batch_size: int = 1, threads: int = 2, force: bool = False) -> dict:
    if batch_size < 1 or threads < 1:
        raise PretrainedError("Batch size and CPU threads must be positive")
    if output.resolve() == fasta.resolve():
        raise PretrainedError("Prediction output cannot overwrite the input FASTA")
    if output.exists() and not force:
        raise PretrainedError("Output already exists; pass --force to replace: " + str(output))
    records = read_fasta(fasta)
    item = manifest["models"][species]
    weights = checkpoint or checkpoint_path(root, item)
    if output.resolve() == weights.resolve():
        raise PretrainedError("Prediction output cannot overwrite a checkpoint")
    verify_checkpoint(weights, item)
    from .strand_model import _lazy_model_imports, build_model, one_hot, spline_basis
    np, torch, _, _ = _lazy_model_imports()
    torch.set_num_threads(threads)
    device = torch.device(device_name)
    if device.type not in ("cpu", "cuda"):
        raise PretrainedError("This prediction workflow supports CPU or CUDA")
    if device.type == "cuda" and not torch.cuda.is_available():
        raise PretrainedError("CUDA is unavailable; use --device cpu")
    # No unsafe-pickle fallback. The published files contain tensors and plain metadata.
    checkpoint_data = torch.load(weights, map_location="cpu", weights_only=True)
    if checkpoint_data.get("species") != SPECIES[species] or checkpoint_data.get("head_hidden") != 256:
        raise PretrainedError("Checkpoint species or head does not match the publication model")
    if (checkpoint_data.get("model_version") != manifest["model_version"] or
            checkpoint_data.get("source_checkpoint_sha256") != item["source_checkpoint_sha256"]):
        raise PretrainedError("Checkpoint model version or source identity differs from its manifest")
    model = build_model(256)
    model.spline.basis = torch.from_numpy(spline_basis(64, 16))
    state = checkpoint_data["state_dict"]
    # Original DataParallel training can leave this fixed, nonlearned buffer empty.
    if state.get("spline.basis") is not None and state["spline.basis"].numel() == 0:
        state = dict(state)
        state["spline.basis"] = model.spline.basis
    model.load_state_dict(state, strict=True)
    del checkpoint_data, state
    model.to(device).eval()
    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.cudnn.benchmark = False
    rows = []
    with torch.inference_mode():
        for start in range(0, len(records), batch_size):
            batch = records[start:start + batch_size]
            tensor = torch.from_numpy(np.stack([one_hot(sequence) for _, sequence in batch])).to(device)
            values = model(tensor).float().cpu().numpy()
            for (identifier, _), pair in zip(batch, values):
                plus, minus = float(pair[0]), float(pair[1])
                if not all(math.isfinite(value) and 0 <= value <= 1 for value in (plus, minus)):
                    raise PretrainedError("Prediction is nonfinite or outside [0,1]: " + identifier)
                rows.append({"sequence_id": identifier, "prediction_plus": plus, "prediction_minus": minus})
    _write_predictions(output, rows)
    return {"species": species, "sequences": len(rows), "device": str(device),
            "precision": "float32_without_autocast", "batch_size": batch_size,
            "checkpoint_sha256": item["sha256"], "fasta_sha256": sha256_file(fasta),
            "torch_version": str(torch.__version__), "output": str(output)}


def read_prediction_table(path: Path) -> dict:
    values = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise PretrainedError("Expected TSV columns: " + ", ".join(FIELDS))
        for row in reader:
            identifier = row["sequence_id"]
            if not identifier or identifier in values:
                raise PretrainedError("Missing or duplicate sequence_id in " + str(path))
            pair = tuple(float(row[field]) for field in FIELDS[1:])
            if not all(math.isfinite(value) and 0 <= value <= 1 for value in pair):
                raise PretrainedError("Expected and predicted scores must be finite and in [0,1]")
            values[identifier] = pair
    if not values:
        raise PretrainedError("Prediction table is empty: " + str(path))
    return values


def compare_predictions(actual_path: Path, expected_path: Path, atol: float = 2e-4,
                        rtol: float = 1e-3) -> dict:
    if not all(math.isfinite(value) and value >= 0 for value in (atol, rtol)):
        raise PretrainedError("Comparison tolerances must be finite and nonnegative")
    actual, expected = read_prediction_table(actual_path), read_prediction_table(expected_path)
    if set(actual) != set(expected):
        raise PretrainedError("Predicted sequence IDs do not exactly match the expected table")
    failures = []
    max_error = 0.0
    max_relative_error = 0.0
    for identifier, reference in expected.items():
        for strand, value, target in zip(("plus", "minus"), actual[identifier], reference):
            error = abs(value - target)
            max_error = max(max_error, error)
            if target:
                max_relative_error = max(max_relative_error, error / abs(target))
            if error > atol + rtol * abs(target):
                failures.append({"sequence_id": identifier, "strand": strand,
                                 "actual": value, "expected": target, "absolute_error": error})
    return {"passed": not failures, "sequences": len(expected), "strand_scores": 2 * len(expected),
            "atol": atol, "rtol": rtol, "max_absolute_error": max_error,
            "max_relative_error_nonzero_reference": max_relative_error, "failures": failures}


def verify_example_files(example_dir: Path, species: str, model_item: dict) -> tuple[Path, Path]:
    with (example_dir / "manifest.json").open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("model_version") != "strand_seiplant_v1" or manifest.get("window_bp") != WINDOW_BP:
        raise PretrainedError("Example manifest has a different model or sequence length")
    item = manifest["species"][species]
    if item["independent_SRX"] != model_item["independent_SRX"]:
        raise PretrainedError("Example cohort does not match the model")
    if item["checkpoint_sha256"] != model_item.get("source_checkpoint_sha256", model_item["sha256"]):
        raise PretrainedError("Example scores were not generated by this source checkpoint")
    paths = []
    for key in ("fasta", "expected"):
        path = (example_dir / item[key]).resolve()
        try:
            path.relative_to(example_dir.resolve())
        except ValueError as exc:
            raise PretrainedError("Example path escapes its directory") from exc
        if sha256_file(path) != item[key + "_sha256"]:
            raise PretrainedError("Example input checksum mismatch: " + str(path))
        paths.append(path)
    records, expected = read_fasta(paths[0]), read_prediction_table(paths[1])
    if len(records) != item["n_sequences"] or {name for name, _ in records} != set(expected):
        raise PretrainedError("Example sequence count or IDs differ from their reference")
    return paths[0], paths[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("download", "predict", "check-example"):
        command = sub.add_parser(name)
        command.add_argument("--root", type=Path, default=repo_root(), help="Repository root")
        command.add_argument("--manifest", type=Path, help="Default: <root>/publication/metadata/pretrained_models.json")
        choices = list(SPECIES) + (["all"] if name != "predict" else [])
        command.add_argument("--species", choices=choices, required=True)
        command.add_argument("--force", action="store_true", help="Replace an existing differing checkpoint or prediction output")
        if name != "download":
            command.add_argument("--device", default="cpu")
            command.add_argument("--batch-size", type=int, default=1)
            command.add_argument("--threads", type=int, default=2)
        if name == "predict":
            command.add_argument("--fasta", type=Path, required=True)
            command.add_argument("--output", type=Path, required=True)
            command.add_argument("--checkpoint", type=Path, help="Alternate location; published hash still required")
        if name == "check-example":
            command.add_argument("--example-dir", type=Path)
            command.add_argument("--output-dir", type=Path)
            command.add_argument("--atol", type=float, default=2e-4,
                                 help="Absolute tolerance for CPU/accelerator numerical differences (default: 0.0002)")
            command.add_argument("--rtol", type=float, default=1e-3,
                                 help="Relative tolerance; pass if error <= atol + rtol*abs(expected) (default: 0.001)")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        root = args.root.resolve()
        manifest = load_manifest(args.manifest or root / "publication/metadata/pretrained_models.json")
        selected = list(SPECIES) if args.species == "all" else [args.species]
        result = {}
        passed = True
        for species in selected:
            if args.command == "download":
                result[species] = download_checkpoint(root, manifest["models"][species], args.force)
            elif args.command == "predict":
                result[species] = predict_fasta(root=root, manifest=manifest, species=species,
                    fasta=args.fasta, output=args.output, checkpoint=args.checkpoint,
                    device_name=args.device, batch_size=args.batch_size, threads=args.threads, force=args.force)
            else:
                example_dir = args.example_dir or root / "examples/real_dna"
                output_dir = args.output_dir or root / "results/pretrained_example"
                fasta, expected = verify_example_files(example_dir, species, manifest["models"][species])
                output = output_dir / (species + ".predictions.tsv")
                if expected.resolve() == output.resolve():
                    raise PretrainedError("Prediction output cannot overwrite the expected table")
                read_prediction_table(expected)  # Fail before loading a large model for a missing/bad reference.
                prediction = predict_fasta(root=root, manifest=manifest, species=species,
                    fasta=fasta, output=output,
                    device_name=args.device, batch_size=args.batch_size, threads=args.threads, force=args.force)
                comparison = compare_predictions(output, expected, args.atol, args.rtol)
                result[species] = {"prediction": prediction, "comparison": comparison}
                passed = passed and comparison["passed"]
                with (output_dir / (species + ".check.json")).open("w", encoding="utf-8") as stream:
                    json.dump(result[species], stream, indent=2)
                    stream.write("\n")
            print(json.dumps({species: result[species]}, indent=2), flush=True)
        return 0 if passed else 1
    except (PretrainedError, OSError, ValueError, KeyError, ImportError, RuntimeError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
