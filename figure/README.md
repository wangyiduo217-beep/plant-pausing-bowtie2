# GROPlant figures and results

This folder keeps the frozen result data, figure previews and plotting scripts
together. The data describe **`strand_seiplant_v1`**, repeats retained,
Arabidopsis / wheat / maize = **14 / 5 / 3 SRX**, covering all **193,508** formal
held-out test windows. Plot-only reproduction does not require trained weights,
FASTQ/BAM files, a full reference genome or a GPU.

## Where to look

| Content | Location |
| --- | --- |
| Ready-to-view result figures | [preview/README.md](preview/README.md) |
| Selected final cross-species radar + tracks, full figure and panels a-f | [plotting/radar_tracks/README.md](plotting/radar_tracks/README.md) |
| Exact inputs for the final radar/track figure | [data/radar_tracks/](data/radar_tracks/) |
| Full window values and region assignments | [data/regional/](data/regional/) |
| Selected genomic tracks | [data/tracks/](data/tracks/) |
| ROC/PR inputs and reference metrics | [data/roc/](data/roc/) |
| CPM-decile estimates, intervals and saved bootstrap draws | [data/deciles/](data/deciles/) |
| Signal-frequency reference tables | [data/frequency/](data/frequency/) |
| Six R figure families and editable plot settings | [plotting/README.md](plotting/README.md) |
| Python ROC/PR, supporting plots and input checks | [scripts/](scripts/) |
| Quantity definitions and statistical conventions | [docs/publication-figures.md](docs/publication-figures.md) |

The R families are complete-test correlations, extended tracks, regional
correlations, dual-axis violins, CPM-decile correlations and signal frequencies.
The selected final cross-species radar/track composite has a separate entry point;
its eight plotting inputs and original A4 reference PDF are now included.
The decile default places both strands side by side in each species panel.
Frequency plots include linear-percent and logarithmic-y views.

## Redraw the results

Run commands from the **repository root**. Install R 4.2.1 or later with Cairo
support; the R plots use standard R packages only. Use Arial for the requested
publication typography, or set `GROSEQ_FONT_FAMILY` to another installed font.

```bash
Rscript --vanilla figure/plotting/run_all.R
```

For the selected final composite (not included in the legacy six-family runner):

```bash
Rscript --vanilla figure/plotting/radar_tracks/plot_radar_tracks_composite.R
Rscript --vanilla figure/plotting/radar_tracks/plot_components.R --component all
python figure/scripts/verify_radar_tracks_inputs.py
```

For the Python ROC/PR figure:

```bash
python -m pip install -r figure/requirements.txt
python figure/scripts/plot_publication_roc_pr.py --font-family "DejaVu Sans"
```

Choose an installed font appropriate for your system. The shared
[installation guide](../docs/acceptance.md) explains fonts, R paths and the
tested environment. Plotting settings and individual commands are in the
[R plotting guide](plotting/README.md).

Regenerated PDFs, PNGs and numerical tables go to **`figure/results/`**. This
local output directory is created by the scripts and excluded from Git;
the original data and previews remain versioned in `data/` and `preview/`.
Rerunning a plot replaces its generated outputs for the same settings.

## Check the data

After installing the Python dependencies, run from the repository root:

```bash
python figure/scripts/verify_publication_inputs.py
python figure/scripts/verify_decile_inputs.py
```

These checks verify file identities and full-window consistency, and reproduce
the decile estimates and frequency counts. Frozen bootstrap draws are supplied
for verification; redrawing the decile figure does not resample them.
The [model folder](../model/README.md) separately contains the weights and
real-DNA inference workflow. The root `check_installation.py` runs both sets of
checks together.

Plotting code uses MIT. Original project-produced tables, bootstrap arrays and
figures use CC BY 4.0 with attribution **GROPlant contributors**; see the
[shared license scope](../docs/licensing.md).
