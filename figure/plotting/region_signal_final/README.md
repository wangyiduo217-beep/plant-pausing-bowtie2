# Final six-panel regional signal figure

Reproduces the selected **`region_signal_six_panels_A4.pdf`** dated 2026-11-11.
The original selected PDFs are preserved under [preview](../../preview/region_signal_final/).

## Reproduce

From the repository root, using Python 3.12 and installed Arial:

```sh
python -m pip install -r figure/plotting/region_signal_final/requirements.txt
python figure/plotting/region_signal_final/build_region_signal.py
```

For a non-Windows machine, pass `--font-dir /path/to/licensed/Arial/fonts` containing
`arial.ttf`, `arialbd.ttf`, and `ariali.ttf`. Font binaries are not redistributed.
`--output-dir` can select another directory. The default is
`figure/results/region_signal_final/`. Paths are resolved from the script, so the
command also works from another working directory.

Outputs: the A4 PDF with its original legend; the 183 × 170 mm PDF, SVG and
600-dpi PNG; six separate panel PDFs in `components/`; and numerical/layout
checks in `validation/`. Individual PDFs are vector clips of the **newly drawn**
figure. Plotting is performed from the supplied numeric summary, after checking
every summary against the full per-window data. No inference, filtering or
bootstrap sampling is run.

| Panel | Species | Output strand |
|---|---|---|
| a | Arabidopsis thaliana | plus |
| b | Arabidopsis thaliana | minus |
| c | Triticum aestivum | plus |
| d | Triticum aestivum | minus |
| e | Zea mays | plus |
| f | Zea mays | minus |

All panels share `build_region_signal.py`; change neither the underlying data nor
the category order to select a panel. Separate exports include the panel's axes,
eight category names and exact window counts. Species/strand labels are also
recorded in `validation/panel_map.json`.

## Data and interpretation

The three original compressed tables are reused byte-for-byte from
[`figure/data/regional`](../../data/regional/), with 45,066 Arabidopsis,
98,656 wheat and 49,786 maize held-out windows. The selected figure's frozen
48-row summary and original category definition are in
[`figure/data/region_signal_final`](../../data/region_signal_final/).
`data_manifest.json` records every input's SHA-256 and repository-relative path.

Bars use the left axis: Pearson r between `log2(1 + mean SRX CPM)` and
`log2(1 + predicted score)` within each region and strand. Orange points/lines
use the independent right axis: arithmetic mean of the **raw experimental CPM**,
including zero-signal windows. The eight regions use the original midpoint
annotation and canonical transcript reference. They are not the later whole-peak
overlap categories. The model is `strand_seiplant_v1`, with repeats retained and
14/5/3 SRX libraries for Arabidopsis/wheat/maize; wheat is the original five-library
version. See [the unchanged figure legend](figure_legend_EN.txt).

`verify_region_inputs.py` checks source hashes and recomputes all 48 rows,
including Pearson, Spearman, mean, sample SD, median and zero counts.

## R alternative and provenance

```sh
Rscript --vanilla figure/plotting/region_signal_final/run_one_page.R
```

This standard-R alternative recalculates the data and writes to
`figure/results/region_signal_final_R/`. It preserves the statistics and design,
but its Cairo page rounding and text placement differ slightly from the chosen
Python PDF. Use the Python entry above for the selected layout. Original Python
and R sources are retained unchanged under `provenance/`; their historical
machine paths are not the portable execution entry points.
