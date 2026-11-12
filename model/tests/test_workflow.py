"""Command contracts and safeguards; no sequencing tools or downloads required."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.workflow import WorkflowError, run_sample


def stages(plan):
    return {stage["stage"]: stage for stage in plan["stages"]}


def option(argv, name):
    return argv[argv.index(name) + 1]


class WorkflowContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="pausing tests ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "not created"
        self.sample = {"run": "SRR_TEST", "sample": "replicate 1", "reads": [str(self.root / "read 1.fastq")]}
        self.settings = {"reference_index": str(self.root / "reference index"), "trimming": {"approved": True}}

    def plan(self, sample=None, settings=None):
        return run_sample(sample or self.sample, settings or self.settings, self.output, dry_run=True)

    def test_dry_run_does_not_read_inputs_or_create_output(self):
        plan = self.plan()
        self.assertEqual(plan["status"], "dry_run")
        self.assertFalse(self.output.exists())
        self.assertFalse(Path(self.sample["reads"][0]).exists())

    def test_spaces_and_cutadapt_expression_are_single_unchanged_arguments(self):
        self.settings["trimming"]["adapter_r1"] = ["A{10};min_overlap=10"]
        plan = stages(self.plan())
        trim = plan["trimming"]["argv"]
        genome = plan["genome_alignment"]["pipeline"][0]
        self.assertEqual(option(trim, "-a"), "A{10};min_overlap=10")
        self.assertIn(self.sample["reads"][0], trim)
        self.assertEqual(option(genome, "-x"), self.settings["reference_index"])
        self.assertEqual(option(genome, "--rg"), "SM:replicate 1")

    def test_single_end_poly_a_and_no_paired_only_flags(self):
        self.settings["trimming"]["poly_a"] = True
        self.settings["rrna_index"] = str(self.root / "rrna index")
        plan = stages(self.plan())
        trim = plan["trimming"]["argv"]
        rrna = plan["rrna_filter"]["argv"]
        genome = plan["genome_alignment"]["pipeline"][0]
        self.assertIn("--poly-a", trim)
        self.assertEqual(option(trim, "-m"), "25")
        self.assertIn("--un-gz", rrna)
        self.assertEqual(option(rrna, "--un-gz"), "SRR_TEST.nonrrna.fastq.gz")
        self.assertIn("-U", genome)
        for args in (rrna, genome):
            self.assertNotIn("-1", args)
            self.assertNotIn("--dovetail", args)
            self.assertNotIn("--no-mixed", args)

    def test_nextflex_two_passes_and_distinct_rrna_genome_pair_rules(self):
        self.sample["reads"].append(str(self.root / "read 2.fastq"))
        self.settings.update(rrna_index=str(self.root / "rrna"), threads={"rrna": 3, "alignment": 5, "sort": 2})
        self.settings["trimming"].update(random_clip=4, adapter_r1=["AGATCGGAAGAGC"], adapter_r2=["AGATCGGAAGAGC"])
        plan = stages(self.plan())
        trim, clip = plan["trimming"]["argv"], plan["random_clip"]["argv"]
        self.assertEqual(option(trim, "-m"), "33")
        self.assertEqual(option(clip, "-m"), "25")
        for flag in ("-u", "-U"):
            self.assertEqual([clip[i + 1] for i, value in enumerate(clip) if value == flag], ["4", "-4"])
        self.assertEqual(option(trim, "-p"), clip[-1])
        rrna = plan["rrna_filter"]["argv"]
        genome, sort = plan["genome_alignment"]["pipeline"]
        self.assertEqual(option(rrna, "-p"), "3")
        self.assertEqual(option(genome, "-p"), "5")
        self.assertEqual(option(sort, "-@"), "2")
        self.assertIn("--un-conc-gz", rrna)
        # Bowtie2's gzip wrapper cannot safely redirect to paths with spaces.
        # Its output argument stays a safe basename; cwd carries the real path.
        self.assertEqual(option(rrna, "--un-conc-gz"), "SRR_TEST.nonrrna_%.fastq.gz")
        self.assertEqual(Path(plan["rrna_filter"]["cwd"]), Path(option(genome, "-1")).parent)
        self.assertNotIn("--dovetail", rrna)
        self.assertNotIn("-X", rrna)
        self.assertIn("--dovetail", genome)
        self.assertEqual(option(genome, "-X"), "1000")
        for args in (rrna, genome):
            self.assertIn("--very-sensitive", args)
            self.assertIn("--no-mixed", args)
            self.assertIn("--no-discordant", args)

    def test_csi_primary_filter_and_no_deduplication_command(self):
        plan = self.plan()
        by_stage = stages(plan)
        for stage in ("bam_index", "high_confidence_index"):
            self.assertIn("-c", by_stage[stage]["argv"])
        high = by_stage["high_confidence"]["argv"]
        self.assertEqual(option(high, "-q"), "20")
        self.assertEqual(option(high, "-F"), "2820")
        self.assertEqual(2820 & 1024, 0, "The duplicate bit must remain unfiltered")
        self.assertTrue(plan["outputs"]["bam_index"].endswith(".bam.csi"))
        self.assertNotIn("markdup", json.dumps(plan))

    def test_sra_plan_uses_split3_and_keeps_input_path_whole(self):
        sample = {"run": "SRR_TEST", "sra": str(self.root / "raw input.sra"), "layout": "PAIRED"}
        conversion = stages(self.plan(sample))["conversion"]["argv"]
        self.assertEqual(conversion[0], "fasterq-dump")
        self.assertIn(sample["sra"], conversion)
        self.assertIn("--split-3", conversion)

    def test_invalid_identifiers_are_rejected_without_directories(self):
        for bad in ("", "..", "../escape", "/absolute", "a/b", "a\\b", "a\n", "x" * 129):
            with self.subTest(run=bad), self.assertRaises(ValueError):
                self.plan({**self.sample, "run": bad})
        self.assertFalse(self.output.exists())

    def test_invalid_read_layouts_are_rejected(self):
        for change in ({"reads": []}, {"reads": "reads.fastq"}, {"reads": ["a", "b", "c"]},
                       {"reads": ["same", "same"]}, {"layout": "PAIRED"},
                       {"sra": "input.sra"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.plan({**self.sample, **change})
        with self.assertRaises(ValueError):
            self.plan({"run": "SRR_TEST", "sra": "input.sra"})

    def test_invalid_thread_values_are_rejected(self):
        for value in (0, -1, 1.5, True, "4", None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.plan(settings={**self.settings, "threads": {"alignment": value}})
        with self.assertRaises(ValueError):
            self.plan(settings={**self.settings, "threads": {"typo": 4}})

    def test_invalid_trimming_is_rejected(self):
        for change in ({"random_clip": -1}, {"random_clip": True}, {"adapter_r1": "AAAA"},
                       {"adapter_r2": ["AAAA"]}, {"poly_a": "true"}, {"approved": 1}):
            cfg = copy.deepcopy(self.settings)
            cfg["trimming"].update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.plan(settings=cfg)

    @unittest.skipUnless(sys.platform.startswith("linux"), "Real execution guards are Linux-specific")
    def test_unapproved_protocol_rejected_before_creating_output(self):
        self.settings["trimming"]["approved"] = False
        with self.assertRaisesRegex(WorkflowError, "not approved"):
            run_sample(self.sample, self.settings, self.output)
        self.assertFalse(self.output.exists())

    @unittest.skipUnless(sys.platform.startswith("linux"), "Real execution guards are Linux-specific")
    def test_partial_output_is_preserved_and_not_restarted(self):
        # Index placeholders are only sufficient to reach this pre-tool guard.
        # Actual index validity/alignment/completion are tested by smoke_real_tools.py.
        Path(self.sample["reads"][0]).write_text("@read\nACGT\n+\nIIII\n")
        for part in ("1", "2", "3", "4", "rev.1", "rev.2"):
            Path(self.settings["reference_index"] + "." + part + ".bt2").write_bytes(b"guard fixture")
        run_dir = self.output / self.sample["run"]
        run_dir.mkdir(parents=True)
        sentinel = run_dir / "partial result.txt"
        sentinel.write_text("preserve this output")
        with self.assertRaisesRegex(WorkflowError, "partial/unrecognized"):
            run_sample(self.sample, self.settings, self.output)
        self.assertEqual(sentinel.read_text(), "preserve this output")
        self.assertFalse((run_dir / "complete.json").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
