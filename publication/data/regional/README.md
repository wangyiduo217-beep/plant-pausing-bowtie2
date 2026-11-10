# Complete test-window inputs

The three `*.regional_plot_data.tsv.gz` files are lossless real plotting inputs from `strand_seiplant_v1`, with repeats retained and 14/5/3 SRX libraries for Arabidopsis/wheat/maize. They contain 45,066/98,656/49,786 rows. The same files support correlations, regional bars and dual-axis violins.

| Column | Meaning |
| --- | --- |
| `window_index` | Stable 1-based row identifier within a species |
| `chrom` | Reference chromosome name |
| `start0`, `end0` | 0-based half-open genomic interval, exactly 1,024 bp |
| `midpoint0` | `(start0 + end0 - 1) / 2` |
| `region` | Frozen canonical coding-transcript regional assignment |
| `mean_cpm_plus`, `mean_cpm_minus` | Equal-weight mean of independently normalized SRX endpoint CPMs in the full window |
| `prediction_plus`, `prediction_minus` | Frozen model's raw 0-1 regression scores |
| `x_plus`, `x_minus` | `log2(1 + mean_cpm)`; the R scripts recompute these from the raw columns |
| `y_plus`, `y_minus` | `log2(1 + prediction)`; the R scripts recompute these from the raw columns |

`expected_statistics_log2p1.csv` and `expected_regional_statistics_log2p1.tsv` contain the frozen complete-test global and regional correlation/OLS checks. They do not replace recomputation. `region_definition.json` documents the eight region categories, transcript rule, overlap priority, matched reference assemblies, public Ensembl annotation URLs and reference checksums. `data_manifest.json` supplies input file sizes and SHA-256 checksums using publication-relative filenames.

These data come from held-out genomic windows, while the experimental libraries also supplied the training labels in other chromosomes. Agreement with this signal is not an independent experimental validation. A model score is not a CPM prediction or a calibrated probability.
