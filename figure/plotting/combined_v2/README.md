# Selected final workflow, model and evaluation composite

This is the execution entry for **GROseq_combined_v2_A4.pdf** (2026-11-12),
with all 14 panels a–n. [Original final PDFs](../../preview/combined_v2/)
and [exact numeric inputs](../../data/combined_v2/) are included.

## Rebuild the numeric panels, then the complete figure

Use Python 3.12 with installed Arial regular, bold, italic and bold italic.

```sh
python -m pip install -r figure/plotting/combined_v2/requirements-source.txt -r figure/plotting/combined_v2/requirements-layout.txt
python figure/plotting/combined_v2/build_source_panels.py
python figure/plotting/combined_v2/build_combined_v2.py --source-figure2 figure/results/combined_v2_source/GROseq_v1_tracks_OLS_correlations_A4.pdf
```

The first command redraws the three regional track groups and six scatter
panels from all supplied numeric observations. The second places those freshly
drawn panels with the original workflow/model artwork and draws the ROC/PR
curves from their numeric vertices. Both scripts accept `--font-dir` and
`--output-dir`; default input and output paths resolve from the script.
For the commands shown, run from the repository root. Font binaries are not
redistributed. No GPU or model inference is required.

The complete final figure is written to `figure/results/combined_v2/`, along
with its legend, 180 × 238 mm artwork, previews, all **14 standalone panel PDFs**
and validation reports. The default output is excluded from version control.

For quick assembly using the preserved original numeric-panel PDF:

```sh
python figure/plotting/combined_v2/build_combined_v2.py
```

That command reproduces the final layout, but does not recalculate the scatter
or tracks. The two-step workflow above provides the full numeric redraw.
Both paths were validated against the selected final PDF.

## Panel mapping and source boundary

| Panels | Content | Reproduction |
|---|---|---|
| a, b | GRO-seq workflow and model architecture | Exact author-supplied vector PDF plus assembly code |
| c, g, k | Arabidopsis, wheat, maize tracks | `build_source_panels.py`; `data/combined_v2/tracks/` |
| d/e, h/i, l/m | Same species' plus/minus scatter | `build_source_panels.py`; `data/combined_v2/test_windows/` |
| f, j, n | Same species' pooled-strand ROC/PR | `build_combined_v2.py`; `data/combined_v2/curves/` |

For a/b the exact vector asset is available, but its complete original creation
program was not found. These schematic panels have no numeric measurement table;
they are **not** claimed to have a new code-only redraw. Original scatter points
are high-resolution raster layers within the PDF; axes, text and lines retain
their vector representation.

[Detailed layout and panel export documentation](README_LAYOUT.md) ·
[Data-driven source panel methods](README-source-panels.md) ·
[Preserved source/curve checksums](provenance/source_manifest.json) ·
[Numerical equality checks](reference_validation/numerical_equivalence.json).
