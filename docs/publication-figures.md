# Selected figure definitions

All primary figures refer to repeat-retained `strand_seiplant_v1`, with 14 Arabidopsis, **5 wheat**, and 3 maize SRX libraries. The complete formal held-out test windows contain 45,066 / 98,656 / 49,786 records. The two prediction columns always retain genomic plus/minus order.

The experimental libraries also supplied the training labels. These figures evaluate held-out genomic windows/chromosomes; they are not independent-library validation or cross-species transfer experiments.

## Experimental signal versus training target

The regression target used for model training is the confidence-weighted experimental consensus overlap within each 1,024-bp window, divided by 1,024 and rounded to 0.1. Background examples have targets (0, 0). See the model guide for interval construction and support requirements.

The x values in the current correlation/violin figures are a different diagnostic measurement: count eligible aligned SE/read-1 **5-prime endpoints** within the same 1,024-bp window, normalize each SRX by its own full-library eligible endpoint count (CPM), then average SRX CPMs equally. Plot `log2(1 + mean CPM)`. This is not RPKM, per-base coverage, a mean of 128-bp display bins, or the mean of separately logged SRX CPMs.

The correlation y values are `log2(1 + predicted score)`. Pearson r is computed on the displayed transformed pairs; Spearman rho uses their ranks. Ordinary least squares is refitted on the actual plotted population. The primary default retains every test window. Optional score thresholds define a conditional subset and must not replace complete-test evaluation or be selected for a higher test-set correlation.

## Current R figures

See [R plotting entry points and settings](../publication/plotting/README.md). The shared regional tables include exact coordinates, region assignments, raw mean CPM, raw prediction scores, and transformed values. They preserve source numeric precision and original row identity. Counts label windows, not biological replicates.

For dual-axis violins, complete distributions determine density and summaries before the final experimental view is clipped to 0-4. The right prediction axis is independently 0-1. Equal heights between colors do not mean equal signal. The arrows above regions report experimental values exceeding the left-axis limit. This display does not discard zero values or tails and does not use min-max normalization.

Track example coordinates are explicitly recorded with their inputs. Wheat and maize use the updated, right-extended 80-kb displays; Arabidopsis uses 40 kb. These loci were selected for experimental peaks and visible agreement. Whole-test correlations/ROC metrics, rather than selected loci, quantify general performance.

## ROC/PR: positive-label definition

For each held-out 1,024-bp window and each strand, a positive binary truth label means **at least one base** intersects an experimental consensus interval on that same strand. Thus `labels_full == (overlap_bp > 0)`. This is not a threshold on the rounded training score, predicted probability, or experimental CPM. The figure pools plus/minus records in row-major order without dropping observations.

The committed NPZ inputs contain `chrom`, `start0`, `end0`, `predicted` (n x 2), `labels_full` (n x 2), and `overlap_bp` (n x 2), plus matching flattened `truth`/`score` arrays. Loading uses `allow_pickle=False`. All windows remain 1,024 bp.

## HiLearn areas and curve display

`scripts/plot_publication_roc_pr.py` calls `ROC_plot` and `PR_curve` from the unchanged, installed **HiLearn 0.2.2** module. The source SHA-256 is checked. The module is loaded directly to avoid importing unrelated optional HiLearn models, while remaining a normal external dependency rather than a vendored copy.

Each call resets the NumPy seed to **42**. When there are more than 2,000 unique prediction scores, HiLearn selects up to 2,000 score thresholds, evaluating every threshold on **all** records. This is threshold subsampling, not sample subsampling. The native PR display adds (recall, precision) endpoints (1, 0) and (0, 1); these are display conventions, not observed sample precision.

The legend reports native HiLearn **AUROC** and **PR AUC**. Native PR AUC is a trapezoidal integral including those display endpoints; it must not be relabeled AP. The script also computes exact all-threshold rank AUROC and step-integral average precision independently and stores them under separate names in `metrics.tsv`. The full-reference per-strand/micro/macro metrics are included with the inputs. Minor native-versus-full differences are expected and documented, not silently substituted.

The ROC chance reference is the grey diagonal. The native HiLearn PR diagonal option is disabled: the meaningful PR chance reference is positive prevalence, not a diagonal line.

## Runtime and fonts

R figures require standard R and Cairo, with Arial available for the selected typography. The Python ROC/PR dependencies are in `requirements-publication.txt`; `--font-family` and `--font-file` allow an installed font choice. Fonts are not redistributed. The actual family is written in the reproduction record. Default ordinary text is 7 pt and panel letters are 8 pt. PDF/SVG exports preserve editable vector elements; inspect the regenerated output at its final publication size.

## CPM deciles and signal frequency

`publication/plotting/deciles/` redraws the supplied group statistics. For each
species and strand separately, CPM=0 forms a separate group and positive CPM
windows are sorted into ten equal-count groups (sizes differ by at most one).
Ties are resolved by chromosome, start0, end0 and original window index; tied
CPM values may span adjacent groups. Zero-group correlations are NA. Both
experiment and prediction are transformed by log2(1+value) before Pearson r.

Intervals are the 2.5/97.5 percentiles of 1,000 chromosome-stratified paired
100-kb block-bootstrap replicates, with seed 20261110. Occupied physical blocks
are sampled with replacement within chromosomes; one draw supplies weights
for both strands and all groups. Original rank groups remain fixed. The
replicate arrays and definitions are included. Rerendering does not resample.
The grouped layout keeps strands separate; Q1 in opposite strands need not
cover the same absolute CPM interval. Bootstrap intervals describe uncertainty
conditional on the fixed model/test population, not biological replicates or
training/library-selection uncertainty. Range restriction affects within-group r.

`publication/plotting/frequency/` independently bins all raw test-window CPMs
after log2(1+CPM), with fixed width 0.25 over 0–12. Each percentage denominator
is every window in that species/strand panel, including zeros. The first bin
contains both zero and low positive values; the exact zero count is annotated
separately. The log-y detail view changes the display axis only. Percentages
sum to 100%; missing regions are not filled with synthetic zeros.
