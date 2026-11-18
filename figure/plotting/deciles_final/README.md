# Final CPM decile figure: separate strand estimates in grouped panels

This module reproduces the user-selected **GROseq_CPM_deciles_grouped_strands_A4.pdf**. It preserves the original Python figure geometry, green/blue strand colors, Arial 7 pt text, 8 pt panel letters, and the unscaled 183 × 115 mm artwork on an exact A4 page. The unchanged selected PDFs are in [`../../preview/deciles_final/reference/`](../../preview/deciles_final/reference/).

The existing [`../deciles/run_decile_correlations.R`](../deciles/run_decile_correlations.R) also draws these exact statistics. Its layout and caption differ slightly from the selected Python original. It remains available as an R implementation; it should not be described as a pixel-identical copy of this selected final layout.

## Reproduce from the repository root

```sh
python -m pip install -r figure/plotting/deciles_final/requirements.txt
python figure/plotting/deciles_final/run_deciles_final.py
```

The default output directory is `figure/results/deciles_final/`. Input paths are resolved relative to the script, so the command also works from a different working directory when an absolute script path is used. Arial regular, bold and italic fonts must be installed; fonts are not redistributed in this repository. If automatic discovery cannot locate Arial, specify a directory containing `arial.ttf`, `arialbd.ttf`, and `ariali.ttf`:

```sh
python figure/plotting/deciles_final/run_deciles_final.py --font-dir "C:/Windows/Fonts" --output-dir output/deciles_final
```

The default `--panel all` creates the full A4 and 183 × 115 mm figures plus three standalone species panels in PDF, SVG and PNG. Each standalone panel is drawn directly from the same statistics using the same drawing function; it is not a screenshot or raster crop. To export just one panel:

```sh
python figure/plotting/deciles_final/run_deciles_final.py --panel a --output-dir output/decile_a
python figure/plotting/deciles_final/run_deciles_final.py --panel b --output-dir output/decile_b
python figure/plotting/deciles_final/run_deciles_final.py --panel c --output-dir output/decile_c
```

| Panel | Species | Data selection | Output stem |
| --- | --- | --- | --- |
| a | Arabidopsis thaliana | `species=arabidopsis`; plus and minus separately | `panel_a_arabidopsis` |
| b | Triticum aestivum | `species=wheat`; plus and minus separately | `panel_b_wheat` |
| c | Zea mays | `species=maize`; plus and minus separately | `panel_c_maize` |

## Actual plotting data and interpretation

The plot reads the already-published [`../../data/deciles/group_statistics.csv`](../../data/deciles/group_statistics.csv): 66 rows, comprising three species × two strands × 11 groups. The 60 nonzero groups supply `pearson_r`, `ci95_low`, `ci95_high`, and `n_windows`; the six zero-CPM groups supply counts and retain undefined correlations. No new statistics, bootstrap samples, or model predictions are computed.

The following distinctions matter:

- Plus and minus strands are **independently ranked and independently correlated**, then displayed side by side. They are not pooled.
- Positive-CPM windows retain their original Q1–Q10 equal-count rank groups within each species and strand. The same Q label does not imply the same absolute CPM cutoff in both strands. Original `cpm_min` and `cpm_max` are in the CSV.
- Pearson correlations use `log2(1 + mean SRX CPM)` and `log2(1 + original predicted score)` within each group. The zero-CPM group has a constant experimental signal, hence r is undefined.
- Whiskers are the previously calculated 95% intervals from 1,000 chromosome-stratified paired 100 kb block-bootstrap replicates. The plotting script only reads those intervals. The original `.block_bootstrap.npz` files and the full window-level input tables are already included in the repository; see [`../../data/deciles/README.md`](../../data/deciles/README.md).
- This is `strand_seiplant_v1`, with repeat sequences retained and Arabidopsis / wheat / maize using 14 / 5 / 3 independent SRX libraries. The wheat result uses the original five-library model.

The three CSVs shipped with the selected final result are byte-identical to the repository's existing copies. Only `group_statistics.csv` supplies plotted values; `overall_correlation_reference.csv` and `zero_signal_prediction_summary.csv` provide context and are copied into generated output for provenance. No duplicate input data are added by this module.

[`panel_map.json`](panel_map.json) records script/data mappings and SHA-256 values. [`provenance/`](provenance/) retains the unchanged original Python and R source files as provenance, including their original environment assumptions. Use `run_deciles_final.py` for the portable selected-layout reproduction.

## Validation scope

The packaged reproduction was checked against the selected final PDF by rendering both with the same Poppler settings and comparing pixels. Plotted r values, intervals and counts were also checked against the original figure's saved artist values. The PDF files contain embedded Arial fonts and vector graphics. PDF byte hashes can differ because of creation metadata even when rendered output is identical. New standalone panels have their own 70 × 105 mm pages; the complete figure retains the selected original layout.
