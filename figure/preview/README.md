# Current figure previews

All previews use repeat-retained `strand_seiplant_v1`, Arabidopsis / wheat / maize = 14 / 5 / 3 SRX. The full editable figures are regenerated under `figure/results/`.

## Complete-test correlations

![Complete-test correlations](correlations.png)

## HiLearn native ROC and PR curves

![HiLearn ROC and PR](roc_pr.png)

Native PR AUC is not AP. See the figure definitions and the separately stored exact metrics.

## Paired dual-axis violins (Arabidopsis example)

![Paired violins](arabidopsis_dual_axis_violins.png)

The entry point produces all three species, with both strands. Green reads against the left axis; orange against the independent right axis. Zooming does not discard high-tail observations.

## CPM deciles with grouped strands

![Strand-specific CPM-decile correlations](deciles_grouped.png)

The plus/minus strands remain separately ranked and calculated. Points and 95%
intervals use the original frozen statistics; this layout does not pool strands.
See `figure/data/deciles/README.md` for the definitions and draw archives.

## Experimental signal frequencies

![Complete experimental signal frequencies](signal_frequency.png)

All formal test windows, including zeros, contribute to each panel's denominator.
The first bin also contains low positive values; exact zeros are annotated.

![Experimental signal frequencies with logarithmic y axis](signal_frequency_logY.png)

This detail view changes only the y-axis display. Counts, percentages and total
window denominators are identical to the linear-percent view.
