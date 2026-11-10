"""Verify file identities and pair complete ROC inputs with the regional figure tables."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hashes-only", action="store_true", help="Standard-library file-identity checks only.")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "publication/metadata/release.json").read_text(encoding="utf-8"))
    assert manifest["repeat_sequences_retained"] is True
    assert manifest["independent_SRX"] == [14, 5, 3]
    for relative, expected in manifest["files_sha256"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Identity mismatch: {relative}")
    if args.hashes_only:
        print(f"Verified {len(manifest['files_sha256'])} selected file identities.")
        return
    import numpy as np
    total = 0
    for key, n in zip(["arabidopsis", "wheat", "maize"], [45066, 98656, 49786]):
        table = ROOT / "publication/data/regional" / f"{key}.regional_plot_data.tsv.gz"
        with gzip.open(table, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == n
        assert [int(r["window_index"]) for r in rows] == list(range(1, n + 1))
        with np.load(ROOT / "publication/data/roc" / f"{key}.full_input.npz", allow_pickle=False) as data:
            assert data["predicted"].shape == data["labels_full"].shape == (n, 2)
            np.testing.assert_array_equal(data["labels_full"], data["overlap_bp"] > 0)
            np.testing.assert_array_equal(data["truth"], data["labels_full"].ravel(order="C"))
            np.testing.assert_array_equal(data["score"], data["predicted"].ravel(order="C").astype(float))
            np.testing.assert_array_equal(data["chrom"], [r["chrom"] for r in rows])
            np.testing.assert_array_equal(data["start0"], [int(r["start0"]) for r in rows])
            np.testing.assert_array_equal(data["end0"] - data["start0"], np.full(n, 1024))
            for index, strand in enumerate(["plus", "minus"]):
                score = np.array([float(r[f"prediction_{strand}"]) for r in rows])
                cpm = np.array([float(r[f"mean_cpm_{strand}"]) for r in rows])
                np.testing.assert_array_equal(score, data["predicted"][:, index].astype(float))
                assert np.all(np.isfinite(cpm)) and np.all(cpm >= 0)
                assert np.all(np.isfinite(score)) and np.all((score >= 0) & (score <= 1))
                np.testing.assert_allclose(np.log2(1 + cpm),
                    [float(r[f"x_{strand}"]) for r in rows], rtol=0, atol=1e-12)
                np.testing.assert_allclose(np.log2(1 + score),
                    [float(r[f"y_{strand}"]) for r in rows], rtol=0, atol=1e-12)
        total += n
    print(f"Verified {total:,} complete windows, 14/5/3 SRX, repeat-retained v1, matching ROC and regional inputs.")


if __name__ == "__main__":
    main()
