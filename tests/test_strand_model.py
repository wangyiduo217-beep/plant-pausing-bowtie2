"""Contracts for DNA-to-strand-signal model preparation."""

from __future__ import annotations

import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from plant_pausing_bowtie2.strand_model import (
    IndexedFasta,
    load_model_config,
    main,
    parse_region,
    prepare_manifests,
)


def write_fasta(path: Path, records: dict[str, str]) -> None:
    offset = 0
    index = []
    with path.open("wb") as stream:
        for name, sequence in records.items():
            header = f">{name}\n".encode()
            stream.write(header); offset += len(header)
            index.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence)+1}\n")
            data = (sequence + "\n").encode()
            stream.write(data); offset += len(data)
    Path(str(path) + ".fai").write_text("".join(index), encoding="utf-8")


class StrandModelPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="strand model test ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config_path = self.root / "configs" / "model.json"
        self.config_path.parent.mkdir()
        reference = self.root / "reference.fa"
        write_fasta(reference, {chrom: ("ATCG" * 2048) for chrom in ("1", "3", "4")})
        labels = self.root / "positive.tsv.gz"
        with gzip.open(labels, "wt", encoding="utf-8", newline="") as stream:
            fields = ["species", "chrom", "start0", "end0", "y_plus", "y_minus",
                      "binary_plus", "binary_minus", "max_support_plus",
                      "max_support_minus", "dataset_status"]
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for chrom, plus, minus in (("1", 0.5, 0.0), ("3", 0.0, 0.7), ("4", 0.2, 0.3)):
                writer.writerow({"species": "Test plant", "chrom": chrom, "start0": 0,
                                 "end0": 1024, "y_plus": plus, "y_minus": minus,
                                 "binary_plus": int(plus > 0), "binary_minus": int(minus > 0),
                                 "max_support_plus": 2, "max_support_minus": 2,
                                 "dataset_status": "test"})
        self.raw = {
            "output": "../results/model",
            "settings": {"negative_ratio": 1.0, "max_ambiguous_fraction": 0.0},
            "training": {"epochs": 2, "batch_size": 1, "num_workers": 0},
            "species": {"Test plant": {
                "slug": "test_plant", "fasta": "../reference.fa",
                "positive_y": "../positive.tsv.gz",
                "chromosomes": {"train": ["1"], "validation": ["4"], "test": ["3"]}
            }}
        }
        self.config_path.write_text(json.dumps(self.raw), encoding="utf-8")

    def test_plan_is_read_only_and_resolves_paths(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(main(["plan", str(self.config_path)]), 0)
        plan = json.loads(output.getvalue())
        self.assertEqual(plan["target_columns"], ["y_plus", "y_minus"])
        self.assertEqual(plan["window_bp"], 1024)
        self.assertFalse((self.root / "results").exists())

    def test_indexed_fasta_fetches_across_line_contract(self):
        config = load_model_config(self.config_path)
        fasta = IndexedFasta(config["species"]["Test plant"]["fasta"])
        try:
            self.assertEqual(fasta.fetch("1", 1, 9), "TCGATCGA")
        finally:
            fasta.close()

    def test_prepare_preserves_chromosome_splits_and_adds_background(self):
        config = load_model_config(self.config_path)
        summary = prepare_manifests(config)["Test plant"]
        for split in ("train", "validation", "test"):
            self.assertEqual(summary["splits"][split]["positive"], 1)
            self.assertEqual(summary["splits"][split]["background"], 1)
            self.assertEqual(summary["splits"][split]["total"], 2)
            manifest = Path(summary["splits"][split]["manifest"])
            cache = Path(summary["splits"][split]["sequence_cache"])
            self.assertEqual(cache.stat().st_size, 2 * 1024)
            with gzip.open(manifest, "rt", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            expected_chrom = {"train": "1", "validation": "4", "test": "3"}[split]
            self.assertEqual({row["chrom"] for row in rows}, {expected_chrom})
            self.assertEqual({row["source"] for row in rows}, {"positive", "background"})

    def test_overlapping_chromosome_splits_are_rejected(self):
        self.raw["species"]["Test plant"]["chromosomes"]["test"] = ["1"]
        self.config_path.write_text(json.dumps(self.raw), encoding="utf-8")
        with self.assertRaises(ValueError):
            load_model_config(self.config_path)

    def test_region_parser_accepts_commas_and_rejects_invalid_coordinates(self):
        self.assertEqual(parse_region("3:9,141,194-9,207,074"), ("3", 9141194, 9207074))
        with self.assertRaises(ValueError):
            parse_region("3:20-10")


if __name__ == "__main__":
    unittest.main(verbosity=2)
