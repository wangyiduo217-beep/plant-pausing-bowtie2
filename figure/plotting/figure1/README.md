# Figure 1: GRO-seq workflow, model and evaluation

This is the execution entry for **GROseq_combined_v2_A4.pdf** (2026-11-12),
with all 14 panels a–n. [Original final PDFs](../../preview/figure1/)
and [exact numeric inputs](../../data/figure1/) are included.

## Rebuild the numeric panels, then the complete figure

Use Python 3.12 with installed Arial regular, bold, italic and bold italic.

```sh
python -m pip install -r figure/plotting/figure1/requirements-source.txt -r figure/plotting/figure1/requirements-layout.txt
python figure/plotting/figure1/build_source_panels.py
python figure/plotting/figure1/build_combined_v2.py --source-figure2 figure/results/figure1_source/GROseq_v1_tracks_OLS_correlations_A4.pdf
```

After dependency installation, the first plotting command redraws the three regional track groups and six scatter
panels from all supplied numeric observations. The assembly command places those freshly
drawn panels with the original workflow/model artwork and draws the ROC/PR
curves from their numeric vertices. Both scripts accept `--font-dir` and
`--output-dir`; default input and output paths resolve from the script.
For the commands shown, run from the repository root. Font binaries are not
redistributed. No GPU or model inference is required.

The complete final figure is written to `figure/results/figure1/`, along
with its legend, 180 × 238 mm artwork, previews, all **14 standalone panel PDFs**
and validation reports. The default output is excluded from version control.

For quick assembly using the preserved original numeric-panel PDF:

```sh
python figure/plotting/figure1/build_combined_v2.py
```

That command reproduces the final layout, but does not recalculate the scatter
or tracks. The two-step workflow above provides the full numeric redraw.
Both paths were validated against the selected final PDF.

After relocation into `figure1/`, the numerical redraw and complete assembly
were executed again using the default input/output paths. The source-panel PDF
and final A4 PDF matched their selected originals pixel for pixel at 160 dpi;
all 14 panel PDFs were produced. The relocated standalone ROC/PR program also
recomputed all six saved curve arrays exactly for each species. See the
[relocation validation](reference_validation/relocation_reproduction_20261118.json).

## Panel mapping and source boundary

| Panels | Content | Reproduction |
|---|---|---|
| a, b | GRO-seq workflow and model architecture | Exact author-supplied vector PDF plus assembly code |
| c, g, k | Arabidopsis, wheat, maize tracks | `build_source_panels.py`; [`figure/data/figure1/tracks/`](../../data/figure1/tracks/) |
| d/e, h/i, l/m | Same species' plus/minus scatter | `build_source_panels.py`; [`figure/data/figure1/test_windows/`](../../data/figure1/test_windows/) |
| f, j, n | Same species' pooled-strand ROC/PR | `build_combined_v2.py`; [`figure/data/figure1/curves/`](../../data/figure1/curves/) |

For a/b the exact vector asset is available, but its complete original creation
program was not found. These schematic panels have no numeric measurement table;
they are **not** claimed to have a new code-only redraw. Original scatter points
are high-resolution raster layers within the PDF; axes, text and lines retain
their vector representation.

[Detailed layout and panel export documentation](README_LAYOUT.md) ·
[Data-driven source panel methods](README-source-panels.md) ·
[ROC/PR reconstruction from labels and scores](roc_pr/README.md) ·
[Preserved source/curve checksums](provenance/source_manifest.json) ·
[Numerical equality checks](reference_validation/numerical_equivalence.json).
