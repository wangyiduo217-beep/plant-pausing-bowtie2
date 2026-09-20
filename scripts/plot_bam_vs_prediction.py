#!/usr/bin/env python3
"""Plot strand-resolved BAM coverage against DNA-model predictions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare mean per-SRX BAM CPM coverage with model predictions"
    )
    parser.add_argument("y_config", help="GRO-seq library configuration JSON")
    parser.add_argument("model_config", help="Strand-model configuration JSON")
    parser.add_argument("--species", required=True, help="Species name or configured slug")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--output", required=True, help="Output prefix")
    return parser.parse_args()


def _prediction_rows(plus_path: Path, minus_path: Path, chrom: str,
                     start: int, end: int) -> list[dict]:
    def read(path: Path) -> dict[tuple[int, int], float]:
        values = {}
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if fields[0] != chrom:
                    continue
                left, right = int(fields[1]), int(fields[2])
                if right > start and left < end:
                    values[(left, right)] = float(fields[3])
        return values

    plus, minus = read(plus_path), read(minus_path)
    if plus.keys() != minus.keys() or not plus:
        raise RuntimeError("Plus/minus prediction tracks are empty or have different coordinates")
    return [
        {"chrom": chrom, "start0": left, "end0": right,
         "predicted_plus": plus[(left, right)],
         "predicted_minus_signed": minus[(left, right)]}
        for left, right in sorted(plus)
    ]


def _library_coverage(library: dict, chrom: str, start: int, end: int):
    import numpy as np
    import pysam

    length = end - start
    plus_diff = np.zeros(length + 1, dtype=np.float64)
    minus_diff = np.zeros(length + 1, dtype=np.float64)
    denominator = 0.0
    regional_reads = 0
    paired = library["layout"] == "PAIRED"
    for bam_name in library["bams"]:
        bam_path = Path(bam_name)
        if not bam_path.is_file():
            raise FileNotFoundError(bam_path)
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # The production BAM contains primary MAPQ>=20 records. For paired
            # libraries the analysis uses R1 only, so mapped/2 is the indexed
            # library-size estimate used for CPM normalization.
            denominator += bam.mapped / 2.0 if paired else bam.mapped
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if paired and not read.is_read1:
                    continue
                regional_reads += 1
                rna_is_reverse = bool(read.is_reverse) ^ bool(library["reverse_strand"])
                difference = minus_diff if rna_is_reverse else plus_diff
                for block_start, block_end in read.get_blocks():
                    left, right = max(block_start, start), min(block_end, end)
                    if right > left:
                        difference[left - start] += 1.0
                        difference[right - start] -= 1.0
    if denominator <= 0:
        raise RuntimeError(f"No indexed mapped reads for {library['srx']}")
    scale = 1_000_000.0 / denominator
    return (np.cumsum(plus_diff[:-1]) * scale,
            np.cumsum(minus_diff[:-1]) * scale,
            denominator, regional_reads)


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from plant_pausing_bowtie2.strand_model import load_model_config, select_species
    from plant_pausing_bowtie2.y_labels import load_y_config

    args = _arguments()
    if args.start < 0 or args.end <= args.start:
        raise ValueError("Require 0 <= --start < --end")
    y_config = load_y_config(args.y_config)
    model_config = load_model_config(args.model_config)
    species_name, model_item = select_species(model_config, [args.species])[0]
    libraries = [item for item in y_config["libraries"] if item["species"] == species_name]
    if not libraries:
        raise RuntimeError(f"No BAM libraries configured for {species_name}")

    root = Path(model_config["output"]) / model_item["slug"]
    rows = _prediction_rows(
        root / "predictions" / f"{model_item['slug']}.plus.bedGraph",
        root / "predictions" / f"{model_item['slug']}.minus_signed.bedGraph",
        args.chrom, args.start, args.end,
    )

    observed_plus = np.zeros(args.end - args.start, dtype=np.float64)
    observed_minus = np.zeros_like(observed_plus)
    provenance = []
    for library in libraries:
        plus, minus, denominator, regional_reads = _library_coverage(
            library, args.chrom, args.start, args.end
        )
        observed_plus += plus
        observed_minus += minus
        provenance.append({
            "srx": library["srx"], "runs": library["runs"],
            "layout": library["layout"], "reverse_strand": library["reverse_strand"],
            "indexed_eligible_read_estimate": denominator,
            "eligible_reads_overlapping_region": regional_reads,
            "bams": library["bams"],
        })
    observed_plus /= len(libraries)
    observed_minus /= len(libraries)

    for row in rows:
        left = max(row["start0"], args.start) - args.start
        right = min(row["end0"], args.end) - args.start
        row["observed_plus_mean_cpm"] = float(observed_plus[left:right].mean())
        row["observed_minus_signed_mean_cpm"] = -float(observed_minus[left:right].mean())

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = output.with_suffix(".tsv")
    with tsv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    x = np.asarray([(row["start0"] + row["end0"]) / 2 for row in rows])
    obs_plus = np.asarray([row["observed_plus_mean_cpm"] for row in rows])
    obs_minus = np.asarray([row["observed_minus_signed_mean_cpm"] for row in rows])
    pred_plus = np.asarray([row["predicted_plus"] for row in rows])
    pred_minus = np.asarray([row["predicted_minus_signed"] for row in rows])
    figure, axes = plt.subplots(2, 1, figsize=(12, 4.8), sharex=True, constrained_layout=True)
    axes[0].fill_between(x, obs_plus, 0, step="mid", color="#087f5b", alpha=0.85,
                         label="plus BAM coverage")
    axes[0].fill_between(x, obs_minus, 0, step="mid", color="#7048e8", alpha=0.85,
                         label="minus BAM coverage")
    axes[0].set_ylabel("Mean coverage\n(CPM per base)")
    axes[0].set_title(f"{species_name}  {args.chrom}:{args.start:,}-{args.end:,}")
    axes[1].fill_between(x, pred_plus, 0, step="mid", color="#20c997", alpha=0.8,
                         label="plus prediction")
    axes[1].fill_between(x, pred_minus, 0, step="mid", color="#9775fa", alpha=0.8,
                         label="minus prediction")
    axes[1].set_ylabel("Predicted y")
    axes[1].set_xlabel(f"Genomic coordinate on {args.chrom} (bp)")
    for axis in axes:
        axis.axhline(0, color="black", linewidth=0.7)
        axis.set_xlim(args.start, args.end)
        axis.legend(frameon=False, ncol=2, loc="upper right")
        axis.spines[["top", "right"]].set_visible(False)
        bound = max(abs(axis.get_ylim()[0]), abs(axis.get_ylim()[1]))
        axis.set_ylim(-bound, bound)
    png_path, pdf_path = output.with_suffix(".png"), output.with_suffix(".pdf")
    figure.savefig(png_path, dpi=300); figure.savefig(pdf_path); plt.close(figure)

    metadata = {
        "species": species_name, "chrom": args.chrom, "start0": args.start, "end0": args.end,
        "signal": "mean of per-SRX strand-resolved aligned-read coverage, CPM per base",
        "aggregation": "SRRs combined within SRX; SRX libraries receive equal weight",
        "paired_end_rule": "R1 only; indexed mapped-record count divided by two for CPM denominator",
        "prediction_segments": len(rows), "libraries": provenance,
        "outputs": {"png": str(png_path), "pdf": str(pdf_path), "table": str(tsv_path)},
    }
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"outputs": metadata["outputs"], "metadata": str(metadata_path),
                      "libraries": len(libraries)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
