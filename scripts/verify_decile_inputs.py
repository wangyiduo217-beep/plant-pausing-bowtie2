"""Verify rank deciles, frozen bootstrap intervals and complete signal frequencies.

Uses repository data only. Recomputes ranks and point estimates, but reads the
existing 1,000 bootstrap draws rather than rerunning the complete bootstrap.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata

SPECS = [("arabidopsis", 45066), ("wheat", 98656), ("maize", 49786)]
STRANDS = ["plus", "minus"]


def table(path, delimiter=","):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label, tolerance=1e-12):
    error = abs(float(actual) - float(expected))
    require(np.isfinite(error) and error <= tolerance, "{}: {} != {}".format(label, actual, expected))
    return error


def check_repository(root, reconstruct=False):
    deciles, frequencies = root / "publication/data/deciles", root / "publication/data/frequency"
    rows = table(deciles / "group_statistics.csv")
    groups = {(r["species"], r["strand"], int(r["group_id"])): r for r in rows}
    require(len(rows) == len(groups) == 66, "Expected 66 unique group rows")
    rows = table(frequencies / "true_signal_frequency_bins.csv")
    bins = {(r["species"], r["strand"], int(r["bin_index"])): r for r in rows}
    require(len(rows) == len(bins) == 288, "Expected 48 frequency bins per panel")
    rows = table(frequencies / "frequency_summary.csv")
    summaries = {(r["species"], r["strand"]): r for r in rows}
    require(len(rows) == len(summaries) == 6, "Expected six frequency summaries")
    report = dict(passed=False, test_windows=0, nonzero_groups=60,
                  valid_bootstrap_replicates_per_group=1000, frequency_panels=6,
                  max_point_error=0.0, max_CI_error=0.0,
                  first_three_bootstrap_draws_reconstructed=reconstruct,
                  max_reconstructed_bootstrap_error=0.0, species=[])
    seeds = np.random.SeedSequence(20261110).spawn(3)
    for species_index, (species, expected_n) in enumerate(SPECS):
        rows = table(root / "publication/data/regional" / (species + ".regional_plot_data.tsv.gz"), "\t")
        n = len(rows)
        require(n == expected_n, species + ": full-test count")
        chrom = np.array([r["chrom"] for r in rows])
        start = np.array([int(r["start0"]) for r in rows])
        end = np.array([int(r["end0"]) for r in rows])
        index = np.array([int(r["window_index"]) for r in rows])
        require(np.all(end-start == 1024), species + ": window length")
        require(np.array_equal(index, np.arange(1, n+1)), species + ": window order")
        block_ids = np.array(["{}:{}".format(c, b) for c, b in zip(chrom, ((start+end)//2)//100000)])
        labels, inverse = np.unique(block_ids, return_inverse=True)
        chromosomes = np.array([s.rsplit(":", 1)[0] for s in labels])
        memberships, xs, ys = [], [], []
        species_report = dict(species=species, n_windows=n, strands={})
        with np.load(deciles / (species+".block_bootstrap.npz"), allow_pickle=False) as stored:
            bootstrap = stored["pearson_r"]
            require(bootstrap.shape == (1000, 2, 11), species + ": bootstrap shape")
            require(np.isnan(bootstrap[:, :, 0]).all(), species + ": zero group must be NA")
            require(np.isfinite(bootstrap[:, :, 1:]).all(), species + ": nonfinite nonzero draws")
            require(int(stored["B"]) == 1000 and int(stored["block_bp"]) == 100000 and int(stored["root_seed"]) == 20261110, species+": bootstrap settings")
            require(np.array_equal(stored["strand_order"], STRANDS), species+": strand order")
            require(np.array_equal(stored["group_names"], ["0"]+["Q"+str(i) for i in range(1,11)]), species+": group order")
            require(np.array_equal(stored["block_labels"], labels) and np.array_equal(stored["block_chromosome"], chromosomes), species+": block membership")
            for si, strand in enumerate(STRANDS):
                c = np.array([float(r["mean_cpm_"+strand]) for r in rows])
                p = np.array([float(r["prediction_"+strand]) for r in rows])
                require(np.isfinite(c).all() and np.all(c >= 0) and np.isfinite(p).all() and np.all((p >= 0) & (p <= 1)), species+": values")
                x, y = np.log2(1+c), np.log2(1+p)
                # Python tuple sort and quotient/remainder partition, independently
                # implemented from the analysis's lexsort/array_split code.
                positive = [i for i in range(n) if c[i] > 0]
                positive.sort(key=lambda i: (c[i], chrom[i], start[i], end[i], index[i]))
                quotient, remainder = divmod(len(positive), 10)
                membership = np.zeros(n, dtype=int)
                offset, counts = 0, []
                for gi in range(1, 11):
                    count = quotient + (gi <= remainder)
                    membership[positive[offset:offset+count]] = gi
                    counts.append(count)
                    offset += count
                require(offset == len(positive) and max(counts)-min(counts) <= 1, species+": group sizes")
                require(np.array_equal(membership == 0, c == 0), species+": zero grouping")
                memberships.append(membership); xs.append(x); ys.append(y)
                for gi in range(11):
                    selected, row = membership == gi, groups[species, strand, gi]
                    require(int(row["n_windows"]) == int(selected.sum()), species+": group n")
                    require(int(row["n_unique_cpm"]) == len(np.unique(c[selected])), species+": distinct CPM")
                    require(int(row["n_occupied_100kb_blocks"]) == len(np.unique(inverse[selected])), species+": occupied blocks")
                    for field, value in [("cpm_min", c[selected].min()), ("cpm_max", c[selected].max()),
                                         ("cpm_median", np.median(c[selected])), ("log2p1_cpm_sd", x[selected].std(ddof=1)),
                                         ("log2p1_prediction_sd", y[selected].std(ddof=1))]:
                        close(value, row[field], species+": "+field)
                    tie = gi > 1 and c[membership == gi-1].max() == c[selected].min()
                    require(row["tie_at_lower_decile_boundary"].lower() == str(bool(tie)).lower(), species+": tie")
                    if gi == 0:
                        require(all(row[k] == "" for k in ["pearson_r","spearman_rho","ci95_low","ci95_high"]), species+": undefined zero statistics")
                        continue
                    r = np.corrcoef(x[selected], y[selected])[0,1]
                    rho = np.corrcoef(rankdata(c[selected]), rankdata(p[selected]))[0,1]
                    report["max_point_error"] = max(report["max_point_error"], close(r,row["pearson_r"],species+": Pearson"), close(rho,row["spearman_rho"],species+": Spearman"))
                    low, high = np.percentile(bootstrap[:,si,gi], [2.5,97.5])
                    report["max_CI_error"] = max(report["max_CI_error"], close(low,row["ci95_low"],species+": CI low",1e-14), close(high,row["ci95_high"],species+": CI high",1e-14))
                    require(int(row["bootstrap_valid"]) == 1000 and int(row["bootstrap_invalid"]) == 0, species+": bootstrap count")
                require(np.all(x < 12), species+": tail exceeds documented range")
                # Explicit interval assignment includes every zero and positive.
                counts_by_bin = np.bincount(np.floor(x/.25).astype(int), minlength=48)
                require(len(counts_by_bin) == 48 and counts_by_bin.sum() == n, species+": omitted windows")
                frequency = counts_by_bin*100.0/n
                close(frequency.sum(), 100.0, species+": frequency sum")
                for bi, count in enumerate(counts_by_bin, 1):
                    row = bins[species,strand,bi]
                    require(int(row["window_count"]) == int(count) and int(row["n_total"]) == n, species+": bin count")
                    close((bi-1)/4, row["left_log2p1_cpm"], species+": bin left")
                    close(bi/4, row["right_log2p1_cpm"], species+": bin right")
                    close(frequency[bi-1], row["frequency_percent"], species+": bin percentage")
                summary, zero = summaries[species,strand], int((c == 0).sum())
                require(int(summary["n_windows"]) == n and int(summary["n_zero"]) == zero and int(summary["n_positive"]) == n-zero, species+": zero count")
                require(int(summary["n_in_first_bin"]) == int(counts_by_bin[0]) and int(summary["n_positive_in_first_bin"]) == int(counts_by_bin[0])-zero, species+": first bin")
                close(summary["zero_percent"],100*zero/n,species+": zero percentage")
                close(summary["first_bin_percent"],frequency[0],species+": first-bin percentage")
                close(summary["max_log2p1_cpm"],x.max(),species+": tail")
                species_report["strands"][strand] = dict(n_zero=zero,nonzero_group_counts=counts,
                    frequency_count_sum=int(counts_by_bin.sum()),frequency_percent_sum=float(frequency.sum()),
                    max_log2p1_cpm=float(x.max()))
            if reconstruct:
                rng = np.random.default_rng(seeds[species_index])
                weights = np.zeros((50,len(labels)),dtype=int)
                for ch in np.unique(chromosomes):
                    ids = np.flatnonzero(chromosomes == ch)
                    weights[:,ids] = rng.multinomial(len(ids),np.full(len(ids),1/len(ids)),size=50)
                for draw in range(3):
                    selected_rows = np.repeat(np.arange(n),weights[draw,inverse])
                    for si in range(2):
                        for gi in range(1,11):
                            ids = selected_rows[memberships[si][selected_rows] == gi]
                            r = np.corrcoef(xs[si][ids],ys[si][ids])[0,1]
                            error = close(r,bootstrap[draw,si,gi],species+": reconstructed paired draw")
                            report["max_reconstructed_bootstrap_error"] = max(report["max_reconstructed_bootstrap_error"],error)
        report["test_windows"] += n
        report["species"].append(species_report)
    report["passed"] = True
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, help="Optional JSON report")
    parser.add_argument("--reconstruct-bootstrap-draws", action="store_true", help="Also reconstruct first three paired draws, not all 1,000")
    args = parser.parse_args()
    report = check_repository(args.root.resolve(), args.reconstruct_bootstrap_draws)
    text = json.dumps(report, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(text.encode("utf-8"))
    print(text, end="")


if __name__ == "__main__":
    main()
