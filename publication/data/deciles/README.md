# Frozen CPM-decile correlations and block-bootstrap draws

These are the original numerical results for repeat-retained `strand_seiplant_v1`,
with Arabidopsis / wheat / maize = 14 / 5 / 3 independent SRX. `group_statistics.csv`
and the three `*.block_bootstrap.npz` files preserve their archived bytes. The
data manifest records sizes and SHA-256 checksums. The full window inputs are
shared with the other figures at `../regional/*.regional_plot_data.tsv.gz`.

## Group definition and quantities

Each species and each strand is handled independently. CPM = 0 is group 0;
positive-CPM windows are sorted by CPM, then chromosome, start0, end0 and original
window_index, and divided into ten rank groups whose sizes differ by at most one.
Tied CPM values may cross adjacent groups. The two strands are not pooled;
identically named Q groups need not span the same absolute CPM range.

Pearson r uses `log2(1 + mean SRX CPM)` and `log2(1 + original model score)`.
The zero group has a constant observed value, so its correlation and confidence
interval are NA. The 60 positive-CPM groups contain defined estimates.

Confidence intervals are the 2.5th and 97.5th percentiles of 1,000 paired 100-kb
block-bootstrap replicates, stratified by chromosome. Blocks are assigned using
window centers; only occupied test-window blocks enter the sampling frame.
Within each chromosome, each draw resamples the original number of occupied
blocks with replacement. The same draw weights apply to both strands and all
groups. Original group memberships are fixed throughout resampling. Root seed
20261110 is split into species-specific streams with NumPy SeedSequence.

## Files

- `group_statistics.csv`: 66 rows (three species, two strands, zero plus ten
  positive groups), including counts, CPM ranges, Pearson/Spearman, 95% intervals
  and valid bootstrap counts.
- `*.block_bootstrap.npz`: `pearson_r` has shape `(1000, 2, 11)`; use the stored
  `strand_order` and `group_names`. Zero-group draws are NA. The archive also
  records occupied block labels/chromosomes, block size, draw count and root seed.
- `zero_signal_prediction_summary.csv`: counts and predicted-score summaries for
  zero-CPM windows.
- `overall_correlation_reference.csv`: whole-test correlations, distinct from
  the within-decile estimates.
- `high_minus_low_correlation_difference.csv`: exploratory paired-draw contrasts
  against Q1. These are not biological-replicate estimates.

## Reproduce the figure

Run `publication/plotting/deciles/run_decile_correlations.R`. It redraws the
frozen estimates and intervals, without recomputing bootstrap. The default
places the plus/minus strands side by side in each species panel. Set
`GROSEQ_DECILE_LAYOUT=sixpanel` or edit `config.R` for the original six-panel layout.

Range restriction affects within-group correlations. Differences across deciles
are descriptive conditional associations, not causal effects of coverage. The
intervals condition on the fixed model and library selection and do not quantify
retraining or independent-experiment uncertainty.
