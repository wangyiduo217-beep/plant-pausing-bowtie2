"""Single-sample Bowtie2 workflow; real execution requires Linux.

Inputs and references are never changed. Failed/partial runs require inspection
and a new output directory; v1 deliberately does not guess whether to overwrite
or resume them. Fingerprints use file path, size and nanosecond modification
time, not full content hashes of potentially hundreds of gigabytes of input.
"""

from __future__ import annotations

import collections
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from typing import Any


VERSION = "0.1.0"
DEFAULT_THREADS = {
    "conversion": 8, "cutadapt": 8, "rrna": 24, "alignment": 16,
    "sort": 4, "qc": 4, "compression": 8,
}


class WorkflowError(RuntimeError):
    """A run failed a preflight, command or completion check."""


def _utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _write_json(path: Path, data: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _absolute(value: Any, name: str) -> str:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ValueError(f"{name} must be a nonempty filesystem path")
    # abspath does not read files; a dry run also works with nonexistent inputs.
    return os.path.abspath(os.fspath(value))


def _normalize(sample: dict, settings: dict) -> tuple[dict, dict]:
    if not isinstance(sample, dict) or not isinstance(settings, dict):
        raise ValueError("sample and settings must be dictionaries")
    if set(settings) - {"reference_index", "rrna_index", "threads", "trimming", "fastqc", "compress_raw"}:
        raise ValueError("settings contains an unsupported key")
    run = sample.get("run", "")
    if not isinstance(run, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", run):
        raise ValueError("run must be a safe 1-128 character identifier")
    sample_id = sample.get("sample", run)
    if not isinstance(sample_id, str) or not sample_id or any(ord(c) < 32 or ord(c) == 127 for c in sample_id):
        raise ValueError("sample must be a nonempty read-group label without control characters")
    reads, sra = sample.get("reads"), sample.get("sra")
    if bool(reads) == bool(sra):
        raise ValueError("provide exactly one of reads or sra")
    normalized_sample = {"run": run, "sample": sample_id}
    if reads:
        if not isinstance(reads, list) or len(reads) not in (1, 2):
            raise ValueError("reads must contain one single-end or two paired-end FASTQ paths")
        paths = [_absolute(p, "reads") for p in reads]
        if len(set(paths)) != len(paths):
            raise ValueError("paired reads must refer to distinct files")
        layout = "PAIRED" if len(paths) == 2 else "SINGLE"
        if sample.get("layout", layout) != layout:
            raise ValueError("layout disagrees with the number of FASTQ inputs")
        normalized_sample.update(reads=paths, layout=layout)
    else:
        layout = sample.get("layout")
        if layout not in ("SINGLE", "PAIRED"):
            raise ValueError("sra input requires layout SINGLE or PAIRED")
        normalized_sample.update(sra=_absolute(sra, "sra"), layout=layout)
    threads = dict(DEFAULT_THREADS)
    supplied_threads = settings.get("threads", {})
    if not isinstance(supplied_threads, dict) or set(supplied_threads) - set(threads):
        raise ValueError("threads contains an unsupported setting")
    threads.update(supplied_threads)
    if any(type(n) is not int or n < 1 for n in threads.values()):
        raise ValueError("every thread count must be a positive integer")
    trim = {"adapter_r1": [], "adapter_r2": [], "poly_a": False,
            "random_clip": 0, "approved": True}
    supplied_trim = settings.get("trimming", {})
    if not isinstance(supplied_trim, dict) or set(supplied_trim) - set(trim):
        raise ValueError("trimming contains an unsupported setting")
    trim.update(supplied_trim)
    for key in ("adapter_r1", "adapter_r2"):
        if not isinstance(trim[key], list) or any(not isinstance(x, str) or not x or any(ord(c) < 32 or ord(c) == 127 for c in x) for x in trim[key]):
            raise ValueError(f"{key} must be a list of nonempty adapter strings")
    if type(trim["random_clip"]) is not int or trim["random_clip"] < 0:
        raise ValueError("random_clip must be a nonnegative integer")
    for key in ("poly_a", "approved"):
        if type(trim[key]) is not bool:
            raise ValueError(f"{key} must be boolean")
    if normalized_sample["layout"] == "SINGLE" and trim["adapter_r2"]:
        raise ValueError("adapter_r2 is not applicable to single-end input")
    normalized_settings = {
        "reference_index": _absolute(settings.get("reference_index"), "reference_index"),
        "rrna_index": _absolute(settings["rrna_index"], "rrna_index") if settings.get("rrna_index") else None,
        "threads": threads, "trimming": trim,
        "fastqc": settings.get("fastqc", True),
        "compress_raw": settings.get("compress_raw", False),
    }
    for key in ("fastqc", "compress_raw"):
        if type(normalized_settings[key]) is not bool:
            raise ValueError(f"{key} must be boolean")
    return normalized_sample, normalized_settings


def _plan(sample: dict, cfg: dict, work: Path) -> dict:
    run, paired = sample["run"], sample["layout"] == "PAIRED"
    t, trim = cfg["threads"], cfg["trimming"]
    data, qc = work / "data", work / "qc"
    count = 2 if paired else 1

    def fastqs(label: str, parent: Path = data) -> list[str]:
        return [str(parent / f"{run}.{label}{'_' + str(i + 1) if paired else ''}.fastq.gz") for i in range(count)]

    raw = sample.get("reads") or [str(work / "raw" / f"{run}.raw{'_' + str(i + 1) if paired else ''}.fastq") for i in range(count)]
    stages: list[dict] = []

    def command(stage: str, args: list[Any], stdout: str | None = None) -> None:
        entry = {"stage": stage, "argv": [str(x) for x in args]}
        if stdout is not None:
            entry["stdout"] = stdout
        stages.append(entry)

    if "sra" in sample:
        command("conversion", ["fasterq-dump", sample["sra"], "--split-3", "-e", t["conversion"],
                               "-t", work / "tmp", "-O", work / "raw"])
        stages[-1]["postprocess"] = "Verify declared layout; reject singleton reads in paired input; normalize generated filenames."
    clean = fastqs("clean")
    intermediate = fastqs("adapter_trim") if trim["random_clip"] else clean
    ca = ["cutadapt", "-j", t["cutadapt"], "-q", "20", "-m", 25 + 2 * trim["random_clip"],
          "--max-n", "0.1", "--json", qc / "cutadapt.json"]
    if trim["poly_a"]:
        ca += ["--poly-a"]
    for adapter in trim["adapter_r1"]:
        ca += ["-a", adapter]
    for adapter in trim["adapter_r2"]:
        ca += ["-A", adapter]
    if trim["adapter_r1"] or trim["adapter_r2"]:
        ca += ["-O", "8", "-e", "0.1"]
    ca += ["-o", intermediate[0]]
    if paired:
        ca += ["-p", intermediate[1]]
    command("trimming", ca + raw)
    if trim["random_clip"]:
        n = trim["random_clip"]
        ca = ["cutadapt", "-j", t["cutadapt"], "-u", n, "-u", -n,
              "-m", "25", "--json", qc / "random_clip.json", "-o", clean[0]]
        if paired:
            ca += ["-U", n, "-U", -n, "-p", clean[1]]
        command("random_clip", ca + intermediate)
    filtered = clean
    if cfg["rrna_index"]:
        filtered = fastqs("nonrrna")
        bt = ["bowtie2", "-p", t["rrna"], "--very-sensitive", "-x", cfg["rrna_index"], "-S", "/dev/null"]
        if paired:
            bt += ["--no-mixed", "--no-discordant", "-1", clean[0], "-2", clean[1],
                   "--un-conc-gz", f"{run}.nonrrna_%.fastq.gz"]
        else:
            bt += ["-U", clean[0], "--un-gz", Path(filtered[0]).name]
        command("rrna_filter", bt)
        # Bowtie2 2.5.5's wrapper interpolates these gzip output paths into a
        # shell redirection. A safe basename plus cwd also works when the user's
        # parent directory contains spaces; passing an argv list alone does not
        # prevent the wrapper's own unquoted-path bug.
        stages[-1]["cwd"] = str(data)
    bam, high = str(data / f"{run}.sorted.bam"), str(data / f"{run}.primary_mapq20.bam")
    bt = ["bowtie2", "-p", t["alignment"], "--very-sensitive", "--end-to-end", "-x", cfg["reference_index"],
          "--rg-id", run, "--rg", "SM:" + sample["sample"]]
    if paired:
        bt += ["-1", filtered[0], "-2", filtered[1], "--no-mixed", "--no-discordant", "-X", "1000", "--dovetail"]
    else:
        bt += ["-U", filtered[0]]
    sort = ["samtools", "sort", "-@", t["sort"], "-m", "1G", "-T", work / "tmp" / "sort", "-o", bam, "-"]
    stages.append({"stage": "genome_alignment", "pipeline": [[str(x) for x in bt], [str(x) for x in sort]]})
    command("bam_validation", ["samtools", "quickcheck", "-v", bam])
    command("bam_index", ["samtools", "index", "-c", "-@", t["qc"], bam])
    command("flagstat", ["samtools", "flagstat", "-@", t["qc"], "-O", "json", bam], str(qc / "flagstat.json"))
    command("idxstats", ["samtools", "idxstats", bam], str(qc / "idxstats.tsv"))
    command("high_confidence", ["samtools", "view", "-@", t["qc"], "-b", "-q", "20", "-F", "2820", "-o", high, bam])
    command("high_confidence_validation", ["samtools", "quickcheck", "-v", high])
    command("high_confidence_index", ["samtools", "index", "-c", "-@", t["qc"], high])
    if cfg["fastqc"]:
        command("fastqc", ["fastqc", "--threads", min(2, t["qc"]), "--outdir", qc] + clean)
    archives = []
    if cfg["compress_raw"]:
        for i, source in enumerate(raw):
            if not source.endswith(".gz"):
                target = str(work / "raw_archive" / f"{run}.raw{'_' + str(i + 1) if paired else ''}.fastq.gz")
                command("compress_raw", ["pigz", "-p", t["compression"], "-c", source], target)
                command("validate_raw_archive", ["pigz", "-t", target])
                archives.append(target)
    return {"stages": stages, "raw_reads": raw,
            "raw_probe": {"reads_per_file": 100000, "output": str(qc / "raw_probe.json"),
                          "timing": "before trimming, after optional SRA conversion"},
            "outputs": {
        "bam": bam, "bam_index": bam + ".csi", "highconf_bam": high,
        "highconf_index": high + ".csi", "clean_reads": clean,
        "alignment_reads": filtered, "qc_directory": str(qc), "raw_archives": archives,
    }}


def _file_identity(path: str | Path) -> dict:
    p = Path(path).resolve(strict=True)
    if not p.is_file():
        raise WorkflowError(f"Expected a regular file: {p}")
    st = p.stat()
    return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _index_identity(prefix: str) -> list[dict]:
    for extension in ("bt2", "bt2l"):
        paths = [prefix + "." + part + "." + extension for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
        if all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in paths):
            return [_file_identity(p) for p in paths]
    raise WorkflowError(f"No complete nonempty Bowtie2 index found for prefix: {prefix}")


def _fingerprint(sample: dict, cfg: dict, versions: dict) -> tuple[str, dict]:
    inputs = [_file_identity(p) for p in sample.get("reads", [sample.get("sra")])]
    if any(item["size"] == 0 for item in inputs):
        raise WorkflowError("Input files must not be empty")
    if len({item["path"] for item in inputs}) != len(inputs):
        raise WorkflowError("Paired FASTQ paths resolve to the same input file")
    details = {"version": VERSION, "sample": sample, "settings": cfg, "inputs": inputs, "tools": versions,
               "reference_index": _index_identity(cfg["reference_index"])}
    if cfg["rrna_index"]:
        details["rrna_index"] = _index_identity(cfg["rrna_index"])
    return hashlib.sha256(json.dumps(details, sort_keys=True).encode()).hexdigest(), details


def _versions(plan: dict) -> dict:
    names = set()
    for stage in plan["stages"]:
        for argv in stage.get("pipeline", [stage.get("argv")]):
            names.add(argv[0])
    versions = {}
    for name in sorted(names):
        executable = shutil.which(name)
        if not executable:
            raise WorkflowError(f"Required executable is missing from PATH: {name}")
        result = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=30, check=False)
        if result.returncode:
            raise WorkflowError(f"Cannot read the version of {name}")
        versions[name] = {"executable": executable, "version_output": (result.stdout + result.stderr).strip()}
    return versions


def _stop(children: list[subprocess.Popen]) -> None:
    # Every launched child owns a fresh session/group, including tool wrappers.
    # Kill the group so grandchildren cannot keep writing after failure.
    previous = {s: signal.signal(s, signal.SIG_IGN) for s in (signal.SIGINT, signal.SIGTERM)}
    try:
        for child in children:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        deadline = time.monotonic() + 5
        while any(child.poll() is None for child in children) and time.monotonic() < deadline:
            time.sleep(0.05)
        for child in children:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def _execute(stage: dict, log: Path) -> None:
    children: list[subprocess.Popen] = []
    output = None
    with log.open("a", encoding="utf-8") as stream:
        commands = stage.get("pipeline", [stage.get("argv")])
        stream.write("\n" + _utc() + " " + " | ".join(shlex.join(c) for c in commands) + "\n")
        if stage.get("cwd"):
            stream.write("Working directory: " + stage["cwd"] + "\n")
        stream.flush()
        try:
            if "pipeline" in stage:
                first = subprocess.Popen(commands[0], stdout=subprocess.PIPE, stderr=stream, start_new_session=True, cwd=stage.get("cwd"))
                children.append(first)
                second = subprocess.Popen(commands[1], stdin=first.stdout, stdout=stream, stderr=stream, start_new_session=True, cwd=stage.get("cwd"))
                children.append(second)
                first.stdout.close()
                sort_rc = second.wait()
                if sort_rc:
                    raise subprocess.CalledProcessError(sort_rc, commands[1])
                align_rc = first.wait()
                if align_rc:
                    raise subprocess.CalledProcessError(align_rc, commands[0])
            else:
                if stage.get("stdout"):
                    output = open(stage["stdout"], "xb")
                child = subprocess.Popen(commands[0], stdout=output or stream, stderr=stream, start_new_session=True, cwd=stage.get("cwd"))
                children.append(child)
                rc = child.wait()
                if rc:
                    raise subprocess.CalledProcessError(rc, commands[0])
        except BaseException:
            _stop(children)
            raise
        finally:
            if output:
                output.close()
            for child in children:
                if child.stdout:
                    child.stdout.close()


def _converted_reads(work: Path, sample: dict, planned: list[str], identity: dict) -> None:
    files = sorted((work / "raw").glob("*.fastq"))
    if sample["layout"] == "PAIRED":
        first = [p for p in files if p.name.endswith("_1.fastq")]
        second = [p for p in files if p.name.endswith("_2.fastq")]
        if len(first) != 1 or len(second) != 1 or len(files) != 2:
            raise WorkflowError("Paired SRA produced missing mates or additional singleton FASTQ; inspect raw files before proceeding")
        files = first + second
    elif len(files) != 1:
        raise WorkflowError("Single-end SRA did not produce exactly one FASTQ")
    if any(p.stat().st_size == 0 for p in files):
        raise WorkflowError("SRA conversion produced an empty FASTQ")
    originals = [p.name for p in files]
    for source, destination in zip(files, map(Path, planned)):
        if source != destination:
            if destination.exists():
                raise WorkflowError(f"Refusing to overwrite conversion output {destination}")
            source.rename(destination)
    _write_json(work / "raw" / "conversion.done.json", {
        "finished_utc": _utc(), "input": identity, "original_names": originals,
        "files": [_file_identity(p) for p in planned],
    })


def _raw_probe(raw: list[str], output: Path) -> None:
    """Inspect a bounded prefix; this is not a full FASTQ validation/checksum."""
    probes = []
    for file in raw:
        lengths, motifs, prefixes = collections.Counter(), collections.Counter(), collections.Counter()
        n, min_quality, max_quality = 0, 127, 0
        opener = gzip.open if file.endswith(".gz") else open
        with opener(file, "rt", encoding="ascii") as handle:
            for _ in range(100000):
                name = handle.readline()
                if not name:
                    break
                sequence = handle.readline().rstrip("\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\r\n")
                if not name.startswith("@") or not plus.startswith("+") or not sequence or len(sequence) != len(quality):
                    raise WorkflowError(f"Malformed FASTQ record {n + 1} in {file}")
                lo, hi = min(map(ord, quality)), max(map(ord, quality))
                if lo < 33 or hi > 126:
                    raise WorkflowError(f"Invalid FASTQ quality character in record {n + 1} of {file}")
                n += 1
                lengths[len(sequence)] += 1
                prefixes[sequence[:12]] += 1
                min_quality, max_quality = min(min_quality, lo), max(max_quality, hi)
                for motif in ("TGGAATTCTCGG", "AGATCGGAAGAG", "GATCGTCGGACT", "AAAAAAAAAA", "TTTTTTTTTT"):
                    if motif in sequence:
                        motifs[motif] += 1
        if n == 0:
            raise WorkflowError(f"No FASTQ records found in {file}")
        probes.append({"file": file, "n": n, "lengths": dict(lengths), "motifs": dict(motifs),
                       "top_prefix": prefixes.most_common(10), "quality_ascii_min": min_quality,
                       "quality_ascii_max": max_quality})
    if len(probes) == 2 and probes[0]["n"] != probes[1]["n"]:
        raise WorkflowError("Paired inputs have different record counts within the inspected prefix")
    _write_json(output, {"reads": probes, "max_records_per_file": 100000,
                         "interpretation": "bounded prefix inspection; not full validation, adapter approval, or quality-encoding determination"})


def _remap(value: Any, source: Path, target: Path) -> Any:
    if isinstance(value, str) and value.startswith(str(source) + os.sep):
        return str(target) + value[len(str(source)):]
    if isinstance(value, list):
        return [_remap(v, source, target) for v in value]
    if isinstance(value, dict):
        return {k: _remap(v, source, target) for k, v in value.items()}
    return value


def _validate_products(plan: dict, cfg: dict) -> None:
    output = plan["outputs"]
    files = [output[key] for key in ("bam", "bam_index", "highconf_bam", "highconf_index")]
    files += output["clean_reads"] + output["alignment_reads"] + output["raw_archives"]
    qc = Path(output["qc_directory"])
    files += [str(qc / name) for name in ("raw_probe.json", "cutadapt.json", "flagstat.json", "idxstats.tsv")]
    if cfg["trimming"]["random_clip"]:
        files.append(str(qc / "random_clip.json"))
    if cfg["fastqc"]:
        for read in output["clean_reads"]:
            base = Path(read).name.removesuffix(".fastq.gz")
            files += [str(qc / (base + "_fastqc" + suffix)) for suffix in (".html", ".zip")]
    for file in files:
        if not Path(file).is_file() or Path(file).stat().st_size == 0:
            raise WorkflowError(f"Expected output is missing or empty: {file}")
    for file in qc.glob("*.json"):
        json.loads(file.read_text(encoding="utf-8"))


def run_sample(sample: dict, settings: dict, output_root: Path, dry_run: bool = False) -> dict:
    """Process one sample; return a command plan, completion record, or matched skip.

    ``sample``: run, sample (optional read-group SM), and either reads=[R1, R2]
    /reads=[SE], or sra=<local file> with layout=PAIRED/SINGLE. ``settings``:
    reference_index, optional rrna_index, threads, trimming, fastqc, compress_raw.
    Thread/trimming defaults are in DEFAULT_THREADS and _normalize. Adapter
    expressions are passed verbatim as single Cutadapt arguments.

    Real execution is Linux-only and must run in the main Python thread. No
    partial-run restart/overwrite is automatic. compress_raw makes verified
    gzip copies and preserves all original inputs, including converted FASTQ.
    """
    sample, cfg = _normalize(sample, settings)
    root = Path(_absolute(output_root, "output_root"))
    directory = root / sample["run"]
    work, artifacts = directory / ".work", directory / "artifacts"
    plan = _plan(sample, cfg, work)
    if dry_run:
        return {"status": "dry_run", "run": sample["run"], "sample": sample,
                "settings": cfg, "output_directory": str(directory), **plan}
    if not sys.platform.startswith("linux"):
        raise WorkflowError("Execution requires Linux; dry_run is available on other platforms")
    if threading.current_thread() is not threading.main_thread():
        raise WorkflowError("Execution must run in the main thread so interruption cleanup is reliable")
    if not cfg["trimming"]["approved"]:
        raise WorkflowError("Trimming protocol is not approved; review adapter/end-processing settings first")
    import fcntl

    root.mkdir(parents=True, exist_ok=True)
    root = root.resolve()
    directory = root / sample["run"]
    if directory.is_symlink():
        raise WorkflowError("Refusing a symlink as the sample output directory")
    directory.mkdir(exist_ok=True)
    if directory.resolve().parent != root:
        raise WorkflowError("Sample output directory escapes output_root")
    # Resolve the root before generating commands and provenance paths.
    work, artifacts = directory / ".work", directory / "artifacts"
    plan = _plan(sample, cfg, work)
    lock_path = directory / ".lock"
    descriptor = os.open(lock_path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, "a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise WorkflowError(f"Sample {sample['run']} is already running") from exc
        completed_path = directory / "complete.json"
        if not completed_path.exists():
            unexpected = [p for p in directory.iterdir() if p.name != ".lock"]
            if unexpected:
                raise WorkflowError("Existing partial/unrecognized output is preserved; inspect it and use a new output_root")
        versions = _versions(plan)
        fingerprint, identity = _fingerprint(sample, cfg, versions)
        if completed_path.exists():
            completed = json.loads(completed_path.read_text(encoding="utf-8"))
            if completed.get("fingerprint") != fingerprint:
                raise WorkflowError("Existing completion has different inputs/settings; use a new output_root")
            for item in completed.get("output_identities", []):
                if _file_identity(item["path"]) != item:
                    raise WorkflowError("Completed output has changed; inspect it and use a new output_root")
            if not completed.get("output_identities"):
                raise WorkflowError("Completion record has no validated output identities")
            return {**completed, "status": "skipped", "reason": "matching fingerprint and output identities"}
        if "sra" in sample and shutil.disk_usage(root).free < identity["inputs"][0]["size"] * 18 + 30 * 1024**3:
            raise WorkflowError("Insufficient free space for conservative SRA conversion allowance (18x SRA + 30 GiB)")
        for name in ("data", "qc", "raw", "raw_archive", "tmp", "logs"):
            (work / name).mkdir(parents=True, exist_ok=False)
        state = {"run": sample["run"], "status": "running", "stage": "preflight", "started_utc": _utc(),
                 "fingerprint": fingerprint, "provenance": identity, "tools": versions, "completed_stages": []}
        _write_json(directory / "status.json", state)

        def interrupted(signum: int, frame: Any) -> None:
            raise WorkflowError(f"Interrupted by signal {signum}")

        previous = {s: signal.signal(s, interrupted) for s in (signal.SIGINT, signal.SIGTERM)}
        try:
            for number, stage in enumerate(plan["stages"], 1):
                if stage["stage"] == "trimming":
                    state.update(stage="raw_probe", updated_utc=_utc())
                    _write_json(directory / "status.json", state)
                    _raw_probe(plan["raw_reads"], work / "qc" / "raw_probe.json")
                    state["completed_stages"].append({"stage": "raw_probe", "finished_utc": _utc(),
                                                      "max_records_per_file": 100000})
                state.update(stage=stage["stage"], updated_utc=_utc())
                _write_json(directory / "status.json", state)
                _execute(stage, work / "logs" / f"{number:02d}_{stage['stage']}.log")
                if stage["stage"] == "conversion":
                    _converted_reads(work, sample, plan["raw_reads"], identity["inputs"][0])
                state["completed_stages"].append({"stage": stage["stage"], "finished_utc": _utc(), "command": stage})
            # Detect changed inputs/references while this run was processing.
            if _fingerprint(sample, cfg, versions)[0] != fingerprint:
                raise WorkflowError("An input or reference changed during processing; results remain partial")
            _validate_products(plan, cfg)
            # A fresh directory and no overwrite option guarantee these are ours.
            work.rename(artifacts)
            receipt_path = artifacts / "raw" / "conversion.done.json"
            if receipt_path.exists():
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                _write_json(receipt_path, _remap(receipt, work, artifacts))
            probe_path = artifacts / "qc" / "raw_probe.json"
            probe = json.loads(probe_path.read_text(encoding="utf-8"))
            _write_json(probe_path, _remap(probe, work, artifacts))
            outputs = _remap(plan["outputs"], work, artifacts)
            output_identities = [_file_identity(p) for p in artifacts.rglob("*") if p.is_file() and "tmp" not in p.relative_to(artifacts).parts]
            state.update(status="complete", stage="complete", finished_utc=_utc(), outputs=outputs,
                         output_identities=output_identities,
                         interpretation={"deduplication": "not performed", "high_confidence": "MAPQ >=20; exclude flags 2820; retain organelles",
                                         "rrna_filter": "concordant pairs only for paired input; not exhaustive" if cfg["rrna_index"] else "not performed",
                                         "model_ready": False, "fingerprint_method": "path, size and mtime_ns; no full input content checksum"})
            _write_json(completed_path, state)
            _write_json(directory / "status.json", state)
            return state
        except BaseException as exc:
            state.update(status="failed", error=str(exc), updated_utc=_utc())
            _write_json(directory / "status.json", state)
            raise
        finally:
            for signum, handler in previous.items():
                signal.signal(signum, handler)
