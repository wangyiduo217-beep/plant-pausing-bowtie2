# Reproduce the publication plots

These scripts reproduce the current **`strand_seiplant_v1`** results with repeat sequences retained. The species are *Arabidopsis thaliana*, *Triticum aestivum* and *Zea mays*, with **14, 5 and 3 SRX libraries**, respectively. Their complete held-out test sets contain **45,066, 98,656 and 49,786 windows**. The wheat results here use the original five libraries, not a three-library retraining. N50 model results are not included.

## Requirements

- R 4.2.1 or later, with Cairo graphics support (`capabilities("cairo")`).
- Only packages distributed with R are used, including `graphics`, `grDevices` and `grid`. No CRAN installation, Python process or network access is needed for these R plots.
- Arial must be installed to reproduce the requested typography. Font binaries are not included. To use another installed font, set `Sys.setenv(GROSEQ_FONT_FAMILY = "Liberation Sans")` before sourcing a script, or edit the module's `config.R`. Changing the font can change text layout; inspect the output again.

The defaults use 7 pt ordinary text, 8 pt panel letters and 600 dpi PNG exports. Both a 183 mm wide figure and an A4 proof are exported. The PDFs retain vector graphics.

## Run all six R figure families

From the repository root:

```sh
Rscript --vanilla publication/plotting/run_all.R
```

Or from R/RStudio:

```r
source("publication/plotting/run_all.R", encoding = "UTF-8")
```

When the working directory is elsewhere, source the **full path** to `run_all.R`. Every input and output path is resolved relative to the scripts, not the working directory. Outputs are written under `publication/generated/`, which is excluded from version control. Existing output files for the same settings are replaced when a plot is rerun.

ROC/PR has a separate Python entry point, `scripts/plot_publication_roc_pr.py`; it is not run by this R entry point.

## Individual figure families

| Directory | Entry point | Output directory | Content |
| --- | --- | --- | --- |
| `correlations/` | `run_plots.R` | `publication/generated/correlations/` | All formal test windows, two strands, density-colored points and an OLS line |
| `tracks/` | `run_tracks.R` | `publication/generated/tracks/` | The current extended track examples alongside the complete-test correlations |
| `regional_bars/` | `run_regional_bars.R` | `publication/generated/regional_bars/` | Region-specific Pearson correlation by strand; Spearman is configurable |
| `dual_axis_violins/` | `run_dual_axis_violins.R` | `publication/generated/dual_axis_violins/` | Experimental and predicted distributions grouped within each region, with independent y-axes |
| `deciles/` | `run_decile_correlations.R` | `publication/generated/deciles/` | Separate plus/minus Pearson estimates and 95% block-bootstrap intervals within nonzero CPM deciles; grouped strands by default |
| `frequency/` | `run_signal_frequency.R` | `publication/generated/frequency/` | Complete experimental-signal frequencies, including zeros, with linear-percent and log-y exports |

For example:

```r
source("publication/plotting/dual_axis_violins/run_dual_axis_violins.R",
       encoding = "UTF-8")
```

Each directory contains a small `config.R`. The track script deliberately reuses the correlation/composite drawing engine; it does not maintain a second copy of that implementation. The regional bar implementation reproduces the current coefficients and sample sizes, with its R layout; it is not intended to be pixel-identical to the earlier Python export.

## What the quantities mean

### Correlation plots

One observation is one **1,024 bp genomic window and one strand**. The experimental count is the number of eligible SE/read1 aligned 5-prime endpoints in that same window. Each SRX is normalized independently:

```text
CPM_library = 1,000,000 * window_endpoint_count / eligible_library_endpoint_count
x = log2(1 + equal-weight mean of the SRX CPMs)
y = log2(1 + predicted window score)
```

The model score is a continuous consensus-derived regression output, not predicted CPM. Pearson and Spearman statistics and the ordinary least-squares fit are recomputed from the displayed x/y observations, separately by strand. Both axes start at zero. Density affects point color only; it does not move observations. The line is a fitted regression line, not the identity line.

The main default is **`MIN_PREDICTION_SCORE <- 0`**, retaining every test window, including experimental zeros. The optional threshold comparison:

```r
source("publication/plotting/correlations/run_threshold_comparison.R",
       encoding = "UTF-8")
```

compares raw score cutoffs 0, 0.001, 0.005 and 0.01. It recalculates retained n, correlations and OLS for each subset. These are exploratory conditional views, not a replacement for complete-test evaluation or a validated biological positivity threshold. No experimental-signal cutoff is applied.

### Extended tracks

The current examples use the following 0-based, half-open display intervals:

| Species | Chromosome | Interval | Length |
| --- | --- | --- | --- |
| Arabidopsis | 3 | `[4,610,000, 4,650,000)` | 40 kb |
| Wheat | 3A | `[292,140,000, 292,220,000)` | 80 kb |
| Maize | 6 | `[74,920,000, 75,000,000)` | 80 kb |

Unlike the experimental x values in the formal correlation plots, the **track experimental signal uses 128 bp display bins**, normalized per library to CPM, averaged and then transformed with `log2(1 + mean CPM)`. Predictions use 1,024 bp DNA inputs with a 128 bp sliding step; each raw 0-1 window score is assigned to its central 128 bp for display. This is not base-resolution prediction. The right-hand correlation panels in the composite still describe the complete formal 1,024 bp test windows, not only the selected track interval.

Gray prediction bins are missing predictions caused by invalid reference input, not score zero; their experimental signal remains visible. The current wheat interval contains 24 invalid prediction bins among 624 display bins. The examples were selected for experimental peaks and agreement, so they illustrate patterns rather than providing an unbiased performance estimate.

### Genomic regions and bars

The shared tables already contain the frozen region assignment for every test window. The assignment uses the window midpoint, `(start0 + end0 - 1) / 2`, with one CDS-containing canonical protein-coding transcript per gene from assembly-matched Ensembl Plants release 63. Canonical transcript selection and overlap priority are documented in `publication/data/regional/region_definition.json`, including public annotation URLs and checksums.

The eight categories are upstream 1-3 kb, upstream <1 kb, 5-prime UTR, CDS, intron, 3-prime UTR, downstream <=1 kb and distal. Flanks follow the transcript's strand. Gene-body features take priority over flanks; overlaps are resolved in the documented order. Both model output strands use the same regional annotation. Distal is relative to this coding-transcript reference and is not synonymous with intergenic relative to every noncoding gene.

Bars are single coefficients from one frozen model per species, with n denoting genomic windows. No biological-replicate dots or uncertainty bars are invented. Fewer than three observations or zero variance produces NA rather than a correlation of zero. The scripts check recalculated current values against the frozen reference tables.

### Dual-axis violins

Each species occupies one page with separate positive- and negative-strand panels. Green experimental violins use the **left axis**, `log2(1 + mean SRX CPM)`, displayed from 0 to 4. Orange predicted violins use the **right axis**, `log2(1 + score)`, from 0 to 1. These are different quantities; equal plotted heights across colors do not mean equal absolute signal.

All observations are used to calculate the histogram-based, reflected Gaussian density and the median/interquartile interval. Only the final experimental viewport is clipped to 0-4; high values are not discarded before density estimation. Arrows report the percentage of experimental windows above 4 within that region. A white dot marks the median and a thick line marks the interquartile interval. No min-max scaling or low-prediction filtering is applied. Full-support shape coordinates and recalculated summaries are exported for checking.

## Data and reproducibility

### CPM deciles and grouped strands

The default decile figure is the latest three-panel layout: each species panel
shows plus-strand points at the left and minus-strand points at the right within
each Q1-Q10 category. Each species/strand was ranked independently; the strands
are not pooled and shared Q labels do not imply equal absolute CPM cutoffs.
The 60 Pearson estimates, confidence intervals and group sizes are unchanged.
Zero-CPM counts are shown separately with correlations marked NA.

The script reads the archived `publication/data/deciles/group_statistics.csv`.
It does not rerun the 1,000 chromosome-stratified, paired 100-kb block-bootstrap
replicates. The draw archives and full statistical definitions are included in
that data directory. Both quantities use `log2(1 + value)` before Pearson.

For the original six-panel layout, in R:

```r
Sys.setenv(GROSEQ_DECILE_LAYOUT = "sixpanel")
source("publication/plotting/deciles/run_decile_correlations.R", encoding = "UTF-8")
Sys.unsetenv("GROSEQ_DECILE_LAYOUT")
```

The default grouped export is 183 x 115 mm, the six-panel export 183 x 183 mm;
each also has a separate A4 proof with its explanatory caption. No page scaling
changes the 7-pt text or 8-pt panel letters.

### Experimental signal frequency

The frequency entry reads the complete shared regional window tables and
independently recomputes 48 fixed bins per species/strand panel. The x value is
`log2(1 + mean SRX CPM)`, the bin width is 0.25, and the range is 0-12. The y value
is each bin's window count divided by all windows in the panel, times 100%.
All windows and experimental zeros are retained. The first bin includes low
positive values as well as zeros; the exact zero count is separately annotated.
The recomputed table is checked against the archived 288-bin reference.

The default also exports `true_signal_frequency_logY_detail`, with unchanged
percentages and denominators and a logarithmic y axis. Empty bins are not drawn
on the log axis. Set `MAKE_LOG_Y <- FALSE` in its `config.R` to skip that view.
Both variants export a 183 x 165 mm figure and an A4 captioned proof, as PDF and
600-dpi PNG.

### Shared data

The three complete regional tables are stored once in `publication/data/regional/` and shared by correlations, regional bars, violins and signal-frequency plots. The small extended-track tables are in `publication/data/tracks/`; frozen decile summaries/draws and frequency references are in `publication/data/deciles/` and `publication/data/frequency/`. Those directories contain schemas, relative filenames and SHA-256 manifests. They contain actual frozen plotting inputs, not simulated or subsampled examples. BAM/FASTQ files, genome FASTA files, model checkpoints, font binaries and machine-specific execution logs are not required to reproduce these plots and are not included here.

The scientific data and plotting methods are preserved from the current result scripts. Migration changes are limited to project-relative input/output paths, a shared input representation, configurable installed fonts, current zero-minimum defaults and the all-figure entry point. Generated settings, fitted coordinates, retained row indices and recomputed statistics accompany the outputs. The release's default R workflow is checked from an external working directory in a fresh R process.
