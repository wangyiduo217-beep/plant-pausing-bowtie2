"""Input integrity and complete-reference comparison for the pretrained workflow."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from plant_pausing_bowtie2.pretrained import (
    PretrainedError, SPECIES, checkpoint_path, compare_predictions,
    download_checkpoint, load_manifest, read_fasta, read_prediction_table, verify_example_files,
)


class PretrainedTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def write(self, name, contents):
        path = self.root / name
        path.write_text(contents, encoding="utf-8")
        return path

    def asset(self, contents=b"verified checkpoint"):
        return {"url": "https://example.org/model.pt", "bytes": len(contents),
                "sha256": hashlib.sha256(contents).hexdigest(),
                "checkpoint_relative_path": "results/model/best.pt"}

    def manifest(self):
        counts = (14, 5, 3)
        return {"schema_version": 1, "model_version": "strand_seiplant_v1",
                "models": {slug: dict(self.asset(), independent_SRX=count)
                           for slug, count in zip(SPECIES, counts)}}

    def test_multiline_lowercase_fasta_preserves_ids_and_forward_sequence(self):
        path = self.write("dna.fa", ">window_7 description\natcg\n" + "atcg" * 255 + "\n")
        self.assertEqual(read_fasta(path), [("window_7", "ATCG" * 256)])

    def test_fasta_rejects_ambiguity_short_duplicate_and_headerless_input(self):
        cases = (">a\n" + "A" * 1023 + "N\n", ">a\nACGT\n", "ACGT\n",
                 ">a\n" + "A" * 1024 + "\n>a\n" + "C" * 1024, ">\n" + "A" * 1024, "")
        for contents in cases:
            with self.subTest(contents=contents[:20]):
                with self.assertRaises(PretrainedError):
                    read_fasta(self.write("invalid.fa", contents))

    def test_download_checks_integrity_and_reuses_verified_file(self):
        body = b"verified checkpoint"
        item = self.asset(body)
        with patch("urllib.request.urlopen", return_value=io.BytesIO(body)) as request:
            result = download_checkpoint(self.root, item)
            self.assertEqual(result["status"], "downloaded_and_verified")
            self.assertEqual(checkpoint_path(self.root, item).read_bytes(), body)
            result = download_checkpoint(self.root, item)
            self.assertEqual(result["status"], "already_verified")
            self.assertEqual(request.call_count, 1)

    def test_bad_download_preserves_existing_file_and_removes_partial(self):
        item = self.asset()
        target = checkpoint_path(self.root, item)
        target.parent.mkdir(parents=True)
        target.write_bytes(b"previous file")
        with patch("urllib.request.urlopen", return_value=io.BytesIO(b"bad download")):
            with self.assertRaises(PretrainedError):
                download_checkpoint(self.root, item, force=True)
        self.assertEqual(target.read_bytes(), b"previous file")
        self.assertEqual(list(target.parent.glob("*.partial")), [])

    def test_existing_unverified_file_requires_force(self):
        item = self.asset()
        target = checkpoint_path(self.root, item)
        target.parent.mkdir(parents=True)
        target.write_bytes(b"wrong")
        with patch("urllib.request.urlopen") as request:
            with self.assertRaises(PretrainedError):
                download_checkpoint(self.root, item)
            request.assert_not_called()

    def test_manifest_rejects_wrong_cohort_insecure_url_and_path_escape(self):
        for key, value in (("independent_SRX", 3), ("url", "http://example.org/model.pt"),
                           ("checkpoint_relative_path", "../model.pt"), ("sha256", "bad")):
            manifest = self.manifest()
            manifest["models"]["triticum_aestivum"][key] = value
            with self.subTest(key=key):
                with self.assertRaises(PretrainedError):
                    load_manifest(self.write("manifest.json", json.dumps(manifest)))
        load_manifest(self.write("manifest.json", json.dumps(self.manifest())))

    def test_resolved_checkpoint_path_cannot_escape_root(self):
        item = dict(self.asset(), checkpoint_relative_path="../outside.pt")
        with self.assertRaises(PretrainedError):
            checkpoint_path(self.root, item)

    def test_comparison_joins_by_identifier_and_uses_all_strands(self):
        header = "sequence_id\tprediction_plus\tprediction_minus\n"
        expected = self.write("expected.tsv", header + "a\t0.1\t0.8\nb\t0.2\t0.7\n")
        actual = self.write("actual.tsv", header + "b\t0.2\t0.7\na\t0.100001\t0.8\n")
        result = compare_predictions(actual, expected)
        self.assertTrue(result["passed"])
        self.assertEqual(result["strand_scores"], 4)
        actual.write_text(header + "a\t0.1\t0.2\nb\t0.2\t0.7\n", encoding="utf-8")
        result = compare_predictions(actual, expected)
        self.assertFalse(result["passed"])
        self.assertEqual(result["failures"][0]["strand"], "minus")

    def test_comparison_cannot_pass_a_subset(self):
        header = "sequence_id\tprediction_plus\tprediction_minus\n"
        expected = self.write("expected.tsv", header + "a\t0.1\t0.8\nb\t0.2\t0.7\n")
        actual = self.write("actual.tsv", header + "a\t0.1\t0.8\n")
        with self.assertRaises(PretrainedError):
            compare_predictions(actual, expected)

    def test_nan_or_duplicate_reference_cannot_pass(self):
        header = "sequence_id\tprediction_plus\tprediction_minus\n"
        for content in ("a\tnan\t0.8\n", "a\t0.1\t0.8\na\t0.1\t0.8\n"):
            with self.assertRaises(PretrainedError):
                read_prediction_table(self.write("bad.tsv", header + content))

    def test_example_hashes_and_source_model_are_required(self):
        slug = "arabidopsis_thaliana"
        fasta = self.write(slug + ".fa", ">a\n" + "ACGT" * 256 + "\n")
        expected = self.write(slug + ".expected.tsv", "sequence_id\tprediction_plus\tprediction_minus\na\t0.1\t0.8\n")
        source_sha = "a" * 64
        item = {"independent_SRX": 14, "checkpoint_sha256": source_sha,
                "fasta": fasta.name, "fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
                "expected": expected.name, "expected_sha256": hashlib.sha256(expected.read_bytes()).hexdigest(),
                "n_sequences": 1}
        self.write("manifest.json", json.dumps({"model_version": "strand_seiplant_v1",
                   "window_bp": 1024, "species": {slug: item}}))
        model_item = {"independent_SRX": 14, "source_checkpoint_sha256": source_sha, "sha256": "b" * 64}
        self.assertEqual(verify_example_files(self.root, slug, model_item), (fasta, expected))
        with self.assertRaises(PretrainedError):
            verify_example_files(self.root, slug, dict(model_item, source_checkpoint_sha256="c" * 64))
        fasta.write_text(">a\n" + "TGCA" * 256 + "\n", encoding="utf-8")
        with self.assertRaises(PretrainedError):
            verify_example_files(self.root, slug, model_item)


if __name__ == "__main__":
    unittest.main()
