"""Small real-tool integration test, entirely isolated from production data.

Usage: python tests/smoke_real_tools.py [--work-dir /absolute/new/audit-directory]
Requires Linux with cutadapt, bowtie2, bowtie2-build, samtools, fastqc and pigz.
No downloads/SRA conversion. Two four-record synthetic inputs exercise paired
NEXTflex end trimming, single-end poly-A, rRNA exclusion, alignment and CSI.
This is software validation, not a validation of the project's 42 real samples.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.workflow import DEFAULT_THREADS, WorkflowError, run_sample
from plant_pausing_bowtie2.cli import build_index


def command(argv):
    return subprocess.run(list(map(str, argv)), check=True, capture_output=True, text=True).stdout.strip()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fastq(path, records):
    with path.open("w", encoding="ascii") as stream:
        for name, sequence in records:
            stream.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_fastq(path):
    records = []
    with gzip.open(path, "rt") as stream:
        while name := stream.readline():
            sequence, plus, quality = stream.readline().strip(), stream.readline(), stream.readline().strip()
            assert name.startswith("@") and plus.startswith("+") and len(sequence) == len(quality)
            records.append((name[1:].strip().split()[0], sequence))
    return records


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expect_rejection(sample, settings, output, phrase):
    try:
        run_sample(sample, settings, output)
    except WorkflowError as exc:
        assert phrase in str(exc), str(exc)
    else:
        raise AssertionError("An unsafe reuse was unexpectedly accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, help="New directory to retain the small audit artifacts")
    args = parser.parse_args()
    if not sys.platform.startswith("linux"):
        parser.error("Real-tool smoke testing requires Linux")
    for executable in ("cutadapt", "bowtie2", "bowtie2-build", "samtools", "fastqc", "pigz"):
        if not shutil.which(executable):
            parser.error(f"Required executable not found: {executable}")
    temporary = None
    if args.work_dir:
        audit = args.work_dir.resolve()
        audit.mkdir(parents=True, exist_ok=False)
    else:
        temporary = tempfile.TemporaryDirectory(prefix="pausing_real_smoke_")
        audit = Path(temporary.name)
    started = time.monotonic()
    rng = random.Random(20260912)
    dna = lambda n: "".join(rng.choices("ACGT", k=n))
    genome, rrna, unrelated = dna(4000), dna(700), dna(700)
    genome_fa, rrna_fa = audit / "tiny genome.fa", audit / "tiny rrna.fa"
    genome_fa.write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
    rrna_fa.write_text(">rRNA_synthetic\n" + rrna + "\n", encoding="ascii")
    genome_index, rrna_index = audit / "genome index", audit / "rrna index"
    for fasta, index in ((genome_fa, genome_index), (rrna_fa, rrna_index)):
        build_args = SimpleNamespace(fasta=fasta, prefix=index, threads=1, dry_run=False)
        built = build_index(build_args)
        assert len(built["outputs"]) == 6 and all(Path(path).stat().st_size for path in built["outputs"])
        shard_hashes = {path: digest(path) for path in built["outputs"]}
        try:
            build_index(build_args)
        except FileExistsError:
            pass
        else:
            raise AssertionError("Existing reference index was unexpectedly overwritten")
        assert shard_hashes == {path: digest(path) for path in built["outputs"]}
    settings = {"reference_index": str(genome_index), "rrna_index": str(rrna_index),
                "threads": {name: 1 for name in DEFAULT_THREADS}, "fastqc": True, "compress_raw": True,
                "trimming": {"approved": True, "adapter_r1": ["AGATCGGAAGAGC"],
                             "adapter_r2": ["AGATCGGAAGAGC"], "random_clip": 4, "poly_a": False}}
    # The two genomic pairs share coordinates deliberately: no deduplication.
    fragments = [("genome_1", genome[300:500]), ("genome_2", genome[300:500]),
                 ("rrna", rrna[200:400]), ("unmapped", unrelated[200:400])]
    expected_clean = [[], []]
    mates = [[], []]
    for name, fragment in fragments:
        sequences = (fragment[:60], reverse_complement(fragment[-60:]))
        for i, sequence in enumerate(sequences):
            expected_clean[i].append((name, sequence))
            mates[i].append((name, "ACGT" + sequence + "TGCA" + "AGATCGGAAGAGC"))
    pe_paths = [audit / "paired read 1.fastq", audit / "paired read 2.fastq"]
    for path, records in zip(pe_paths, mates):
        fastq(path, records)
    input_hashes = {str(path): digest(path) for path in pe_paths}
    pe_sample = {"run": "SYNTHETIC_PE", "sample": "synthetic paired replicate", "reads": list(map(str, pe_paths))}
    output = audit / "workflow results"
    pe = run_sample(pe_sample, settings, output)
    assert pe["status"] == "complete"
    for i, path in enumerate(pe["outputs"]["clean_reads"]):
        assert read_fastq(path) == expected_clean[i], "NEXTflex clipping changed insert sequence"
    results = []

    def validate(record, expected_total, expected_high, expected_filtered):
        outputs = record["outputs"]
        command(["samtools", "quickcheck", "-v", outputs["bam"], outputs["highconf_bam"]])
        assert all(Path(outputs[key]).is_file() for key in ("bam_index", "highconf_index"))
        assert not Path(outputs["bam"] + ".bai").exists()
        total = int(command(["samtools", "view", "-c", outputs["bam"]]))
        high = int(command(["samtools", "view", "-c", outputs["highconf_bam"]]))
        assert total == expected_total, (total, expected_total)
        assert high == expected_high, (high, expected_high)
        for path in outputs["alignment_reads"]:
            records = read_fastq(path)
            assert len(records) == expected_filtered
            assert {name for name, _ in records} == {"genome_1", "genome_2", "unmapped"}
        assert record["interpretation"]["model_ready"] is False
        assert record["interpretation"]["deduplication"] == "not performed"
        header = command(["samtools", "view", "-H", outputs["bam"]])
        assert "SO:coordinate" in header and "@RG" in header
        # An indexed region query checks that the CSI is usable, not just present.
        assert int(command(["samtools", "view", "-c", outputs["highconf_bam"], "chrSynthetic:1-4000"])) == high
        assert list(Path(outputs["qc_directory"]).glob("*_fastqc.html"))
        results.append({"run": record["run"], "sorted_records": total, "primary_mapq20_records": high,
                        "rRNA_filtered_records_per_mate": expected_filtered, "CSI_region_query": "passed"})

    validate(pe, expected_total=6, expected_high=4, expected_filtered=3)
    for path in pe_paths:
        assert digest(path) == input_hashes[str(path)], "Original input was modified"
    for original, archive in zip(pe_paths, pe["outputs"]["raw_archives"]):
        with gzip.open(archive, "rb") as stream:
            assert stream.read() == original.read_bytes()
    bam_hash = digest(pe["outputs"]["bam"])
    completed_mtime = (output / pe_sample["run"] / "complete.json").stat().st_mtime_ns
    again = run_sample(pe_sample, settings, output)
    assert again["status"] == "skipped"
    assert digest(pe["outputs"]["bam"]) == bam_hash
    assert (output / pe_sample["run"] / "complete.json").stat().st_mtime_ns == completed_mtime
    changed = copy.deepcopy(settings)
    changed["trimming"]["random_clip"] = 0
    expect_rejection(pe_sample, changed, output, "different inputs/settings")
    # Modify one synthetic output to ensure the completion receipt is checked.
    probe = Path(pe["outputs"]["qc_directory"]) / "idxstats.tsv"
    stat = probe.stat()
    import os
    os.utime(probe, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    expect_rejection(pe_sample, settings, output, "Completed output has changed")
    os.utime(probe, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert run_sample(pe_sample, settings, output)["status"] == "skipped"

    se_path = audit / "single read.fastq"
    # Six non-A bases before the artificial tail avoid an ambiguous A-rich boundary.
    sequences = {}
    for name, sequence in (("genome_1", genome), ("genome_2", genome), ("rrna", rrna), ("unmapped", unrelated)):
        start = next(i for i in range(100, 300) if "A" not in sequence[i + 54:i + 60])
        sequences[name] = sequence[start:start + 60]
    fastq(se_path, [(name, sequence + "A" * 15) for name, sequence in sequences.items()])
    se_settings = copy.deepcopy(settings)
    se_settings["trimming"] = {"approved": True, "poly_a": True, "random_clip": 0}
    se_settings["compress_raw"] = False
    se_sample = {"run": "SYNTHETIC_SE", "sample": "synthetic single replicate", "reads": [str(se_path)]}
    se = run_sample(se_sample, se_settings, output)
    assert read_fastq(se["outputs"]["clean_reads"][0]) == list(sequences.items())
    validate(se, expected_total=3, expected_high=2, expected_filtered=3)
    assert not se["outputs"]["raw_archives"]
    bad_fastq = audit / "malformed input.fastq"
    bad_fastq.write_text("@broken\nACGT\n+\nII\n", encoding="ascii")
    bad_sample = {"run": "SYNTHETIC_BAD", "reads": [str(bad_fastq)]}
    expect_rejection(bad_sample, se_settings, output, "Malformed FASTQ")
    bad_directory = output / bad_sample["run"]
    assert not (bad_directory / "complete.json").exists()
    assert json.loads((bad_directory / "status.json").read_text())["status"] == "failed"
    expect_rejection(bad_sample, se_settings, output, "partial/unrecognized")
    summary = {"status": "passed", "elapsed_seconds": round(time.monotonic() - started, 2), "samples": results,
               "additional_checks": ["CLI reference building", "existing index overwrite rejected", "exact NEXTflex cleaned inserts", "exact poly-A cleaned inserts", "rRNA decoy exclusion",
                                     "coordinate duplicates retained", "original inputs unchanged", "gzip archive roundtrip",
                                     "matching completion skipped", "changed settings rejected", "changed output rejected",
                                     "malformed FASTQ fails without completion", "failed partial result preserved"],
               "limitations": "Tiny synthetic FASTQ integration only; no SRA conversion, production references, scheduler or 42-sample biological validation.",
               "audit_directory": str(audit), "retained": bool(args.work_dir)}
    (audit / "smoke_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if temporary:
        temporary.cleanup()


if __name__ == "__main__":
    main()
