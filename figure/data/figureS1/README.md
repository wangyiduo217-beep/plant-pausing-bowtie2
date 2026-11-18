# Figure S1: exact data for the regional signal figure

The original final figure summary contains 48 rows: 3 species × 2 strands ×
8 regions. `regional_correlation_and_mean_signal.tsv` is the actual plot input.
`previous_regional_correlations.tsv` and `association_across_eight_categories.tsv`
are retained provenance tables and are not new tests added to this figure.
`region_definition.json` documents the original window midpoint annotations.

Full per-window observations are reused from `../regional/*.regional_plot_data.tsv.gz`;
their SHA-256 hashes are identical to those supplied with the selected final PDF.
`data_manifest.json` lists both reused and new inputs. The plotting entry verifies
them and recalculates the summaries before rendering. The 193,508 windows are
counted once per species; both strand columns are analyzed separately.

Main summary fields:

- `species`, `strand`, `region`: group identity.
- `n_windows`: all windows in the category, including zeros.
- `pearson_r`, `spearman_rho`: correlations of the two log2(1 + value) arrays.
- `mean_cpm`, `sd_cpm`, `median_cpm`: raw experimental mean-SRX CPM summaries;
  SD uses `ddof=1`. Only `mean_cpm` is plotted on the right axis.
- `zero_count`, `zero_fraction`: exact experimental zero values.
- `mean_log2p1_cpm`, `log2p1_mean_cpm`: distinct transformation orderings,
  retained for audit; neither replaces the plotted raw CPM mean.

[Code, panel map and commands](../../plotting/figureS1/README.md).
