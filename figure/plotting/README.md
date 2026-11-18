# Reproduce Figures 1, 2, S1 and S2

## Selected final figures: script index

Each top-level subdirectory is one final manuscript figure.
`figure/scripts/` contains data checks, not drawing entry points.

| Figure | Content | Entry files | Data |
|---|---|---|---|
| **Figure 1** | Workflow/model and evaluation, a–n | [figure1/](figure1/README.md): `build_source_panels.py`, then `build_combined_v2.py` | [figure1](../data/figure1/) |
| **Figure 2** | Cross-species radar and tracks, a–f | [figure2/](figure2/README.md): `plot_radar_tracks_composite.R`, `plot_components.R` | [figure2](../data/figure2/) |
| **Figure S1** | Regional correlation and mean signal, a–f | [figureS1/](figureS1/README.md): `build_region_signal.py` | [summary](../data/figureS1/) + [shared windows](../data/regional/) |
| **Figure S2** | CPM deciles with grouped strands, a–c | [figureS2/](figureS2/README.md): `run_deciles_final.py` | [figureS2](../data/figureS2/) |

[View the four complete figures and PDFs](../preview/README.md).

## Environment

Use a separate **Python 3.12 plotting environment**, with the module requirements
below and locally installed Arial regular/bold/italic/bold-italic fonts.
Figure 2 uses **R 4.2.1 or later with Cairo**, with standard R packages only.
Font binaries are not included.

```sh
python -m pip install -r figure/plotting/figure1/requirements-source.txt -r figure/plotting/figure1/requirements-layout.txt -r figure/plotting/figureS1/requirements.txt -r figure/plotting/figureS2/requirements.txt
```

The older model acceptance environment has a different dependency lock; do not
apply it to these final-layout plotting requirements. See the
[installation guide](../../docs/acceptance.md).

## Reproduce all four figures

From the repository root:

```sh
python figure/plotting/reproduce_all.py --font-dir "C:/Windows/Fonts"
```

Use your installed Arial directory on other platforms. Supply `--rscript`
with the full executable path if needed. This entry only redraws figures from
the supplied inputs; it does not download weights or run a model.

## Individual figures and panels

Run after installing the dependencies above:

```sh
# Figure 1: redraw tracks/scatter, then assemble and export all a–n panels
python figure/plotting/figure1/build_source_panels.py
python figure/plotting/figure1/build_combined_v2.py --source-figure2 figure/results/figure1_source/GROseq_v1_tracks_OLS_correlations_A4.pdf

# Figure 2: complete figure, then all separate panels
Rscript --vanilla figure/plotting/figure2/plot_radar_tracks_composite.R
Rscript --vanilla figure/plotting/figure2/plot_components.R --component all

# Figure S1: complete figure and all separate panels
python figure/plotting/figureS1/build_region_signal.py

# Figure S2: complete figure and all separate panels
python figure/plotting/figureS2/run_deciles_final.py
```

For one Figure 2 panel, use `--component a` (or b–f) with
`plot_components.R`. For one Figure S2 panel, add `--panel a` (or b/c).
Figure 1 and S1 export all separate panels by default and have no single-panel
selection flag.

| Figure | Panel mapping | Default regenerated outputs |
|---|---|---|
| 1 | a/b: workflow/model artwork; c/g/k: tracks; d/e, h/i, l/m: scatter; f/j/n: ROC/PR | `figure/results/figure1/`; source redraw in `figure1_source/` |
| 2 | a/b/c: Arabidopsis/maize/wheat radar; d/e/f: their tracks | `figure/results/figure2/figures/` and `components/` |
| S1 | a/b, c/d, e/f: plus/minus for Arabidopsis, wheat, maize | `figure/results/figureS1/` and `components/` |
| S2 | a/b/c: Arabidopsis/wheat/maize, strands estimated separately | `figure/results/figureS2/` |

All output directories are excluded from Git. The original PDF filenames remain
unchanged; the directory identifies the manuscript figure number.

## Supporting sources and checks

Figure 1 a/b reuse the exact vector asset under `figure1/assets/`; the original
primitive-creation program was not recovered. Numeric panels have data-driven
redraw programs. Archived `provenance/original_*` files preserve historical
source bytes and may contain their original paths; use the current entry files
listed above.

[Figure 1 ROC/PR reconstruction](figure1/roc_pr/README.md) retains the original
HiLearn 0.2.2 calculations from complete labels/scores. It is nested under Figure 1
because the same results contribute to that figure.

Shared regional window data remain once under [data/regional/](../data/regional/).
Figure S1 and S2 use these for statistical verification; they are not a separate
result family. The [figure methods](../docs/publication-figures.md) define signal
transforms, labels, regions, metrics and bootstrap intervals.
