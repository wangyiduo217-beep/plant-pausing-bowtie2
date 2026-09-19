"""Contracts for the portable GRO-seq interval and y-label scripts."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from plant_pausing_bowtie2.y_labels import (
    confidence_for_support,
    load_y_config,
    main,
    plan_y_labels,
    rounded_y,
    score_windows,
)


class YLabelTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="y label test ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config_path = self.root / "configs" / "labels.json"
        self.config_path.parent.mkdir()
        self.raw = {
            "output": "../results/y",
            "settings": {"window_bp": 1024, "step_bp": 512, "minimum_support": 2},
            "species": {
                "Test plant": {
                    "slug": "test_plant",
                    "fai": "../references/genome.fa.fai",
                    "mask": "../references/excluded.bed",
                    "dataset_status": "article_analog_primary",
                }
            },
            "libraries": [
                {
                    "srx": "SRX1", "species": "Test plant", "project": "GSE1",
                    "group": "WT", "layout": "PAIRED", "runs": ["SRR1"],
                    "bams": ["../bams/SRR1.primary_mapq20.bam"], "reverse_strand": True,
                },
                {
                    "srx": "SRX2", "species": "Test plant", "project": "GSE1",
                    "group": "WT", "layout": "SINGLE", "runs": ["SRR2"],
                    "bams": ["../bams/SRR2.primary_mapq20.bam"], "reverse_strand": False,
                },
            ],
        }
        self.write()

    def write(self):
        self.config_path.write_text(json.dumps(self.raw), encoding="utf-8")

    def test_config_paths_and_plan_preserve_scientific_parameters(self):
        config = load_y_config(self.config_path)
        self.assertEqual(config["output"], str((self.root / "results/y").resolve()))
        self.assertEqual(config["species"]["Test plant"]["fai"],
                         str((self.root / "references/genome.fa.fai").resolve()))
        plan = plan_y_labels(config)
        paired, single = plan["libraries"]
        paired_view = paired["analysis_inputs"][0]["command"]["pipeline"][0]
        single_view = single["analysis_inputs"][0]["command"]["pipeline"][0]
        self.assertEqual(paired_view[paired_view.index("-f") + 1], "64")
        self.assertNotIn("-f", single_view)
        self.assertIn("-rev", paired["find_peaks"])
        self.assertNotIn("-rev", single["find_peaks"])
        self.assertIn("groseq", paired["find_peaks"])
        self.assertEqual(plan["window_bp"], 1024)
        self.assertEqual(plan["step_bp"], 512)

    def test_plan_cli_does_not_create_output_or_read_bams(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.assertEqual(main(["plan", str(self.config_path)]), 0)
        result = json.loads(stdout.getvalue())
        self.assertEqual(len(result["libraries"]), 2)
        self.assertFalse((self.root / "results").exists())

    def test_invalid_library_configuration_is_rejected(self):
        self.raw["libraries"][1]["srx"] = "SRX1"
        self.write()
        with self.assertRaises(ValueError):
            load_y_config(self.config_path)
        self.raw["libraries"][1]["srx"] = "SRX2"
        self.raw["libraries"][1]["layout"] = "UNKNOWN"
        self.write()
        with self.assertRaises(ValueError):
            load_y_config(self.config_path)

    def test_confidence_and_rounding_match_production_rules(self):
        observed = [2, 3, 4, 4]
        self.assertAlmostEqual(confidence_for_support(2, observed, 14), 0.1)
        self.assertAlmostEqual(confidence_for_support(3, observed, 14), 0.55)
        self.assertAlmostEqual(confidence_for_support(4, observed, 14), 1.0)
        self.assertEqual(confidence_for_support(2, [2, 2], 14), 1.0)
        self.assertEqual(rounded_y(0.049), 0.0)
        self.assertEqual(rounded_y(0.05), 0.1)
        self.assertEqual(rounded_y(0.26), 0.3)
        self.assertEqual(rounded_y(2.0), 1.0)

    def test_strands_are_scored_independently_on_1024_512_grid(self):
        intervals = [
            {"chrom": "chr1", "start": 0, "end": 1024, "strand": "+",
             "confidence": 0.5, "support": 2},
            {"chrom": "chr1", "start": 512, "end": 1536, "strand": "-",
             "confidence": 1.0, "support": 3},
        ]
        scores = score_windows(intervals, {"chr1": 2048}, 1024, 512)
        self.assertEqual(scores[("chr1", 0, 1024)], [0.5, 0.5, 2, 3])
        self.assertEqual(scores[("chr1", 512, 1536)], [0.25, 1.0, 2, 3])
        self.assertEqual(scores[("chr1", 1024, 2048)], [0.0, 0.5, 0, 3])


if __name__ == "__main__":
    unittest.main(verbosity=2)
