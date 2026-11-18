# Final workflow, model and evaluation figure

This package reproduces the selected `GROseq_combined_v2_A4.pdf` dated 2026-11-12. It is the figure with workflow/model diagrams, representative tracks, strand-specific scatter plots, and micro ROC/PR curves. It is distinct from the later cross-species radar-and-tracks figure.

## Run

From the repository root, install the plotting dependencies and run:

```sh
python -m pip install -r figure/plotting/combined_v2/requirements-layout.txt
python figure/plotting/combined_v2/build_combined_v2.py --font-dir /path/to/arial-fonts --output-dir figure/results/combined_v2
```

On Windows, installed Arial fonts default to `C:/Windows/Fonts`; `--font-dir` can be omitted. Other systems need a directory containing `arial.ttf`, `arialbd.ttf`, and `ariali.ttf`. `GROSEQ_FONT_DIR` also sets this directory. Font binaries are not included. The program does not silently substitute a different font.

Optional flags:

- `--source-figure2 /path/to/GROseq_v1_tracks_OLS_correlations_A4.pdf`: assemble using the zero-origin Figure2 rebuilt from the numerical data rather than the preserved original PDF. The input must retain the original A4 physical layout and coordinates.
- `--poppler-preview`: additionally use `pdftoppm` from PATH to render a preview. A PyMuPDF preview is generated without this flag.
- `--no-panels`: skip individual panel exports.

The default output contains the A4 PDF, the 180 x 238 mm artwork PDF, PNG previews, the figure legend, 14 individual panel PDFs, and JSON validation records. The A4 PDF retains the original caption, panel coordinates, data, colors, labels, and scales. Individual panel exports add missing shared labels where needed; they are not used to alter the A4 layout.

## Panel-to-source map

| Panels | Content | Source |
| --- | --- | --- |
| a | GRO-seq workflow | `assets/source_figure1_workflow_model.pdf` |
| b | Sequence-model architecture | Same original vector PDF, plus the original assembler's module outlines and connector adjustments |
| c, g, k | Arabidopsis, wheat, maize tracks | Original zero-origin Figure2 PDF, or the numerical rebuild supplied with `--source-figure2` |
| d, e | Arabidopsis plus/minus scatter | Same Figure2 source |
| h, i | Wheat plus/minus scatter | Same Figure2 source |
| l, m | Maize plus/minus scatter | Same Figure2 source |
| f, j, n | Arabidopsis, wheat, maize micro ROC/PR | `figure/data/combined_v2/curves/{species}.HiLearn_curves.npz` |

Each track panel contains all four displayed tracks. Each scatter export contains only one output strand. Species and sample size are added above standalone scatter panels, and the shared x-axis label is repeated below them. For standalone minus-strand scatter panels, the y-axis label shared in the composite is also repeated. These additions do not affect the final assembled figure.

## What is, and is not, regenerated

The default assembly preserves PDF objects for the workflow/model diagrams and source Figure2 panels, and draws ROC/PR curves directly from saved numerical vertices. It does not rerun training, model inference, read alignment, or peak calling.

Panels a/b are author-supplied vector artwork. The exact original PDF is included, but its complete original creation program was not found. This is an explicit source-asset boundary: the assembly code recreates their placement, annotations, boxes and connectors, rather than claiming to redraw every schematic primitive from a numerical table. An earlier related SVG is deliberately not substituted for the exact final artwork.

The tracks/scatter plotting program and numerical data are packaged separately alongside this assembler. Supplying that program's rebuilt zero-origin A4 PDF with `--source-figure2` uses the new numerical redraw. The original scatter points were rasterized at high resolution inside a PDF; axes, labels and other vector elements remain editable. The figure is therefore a hybrid PDF, not a screenshot redraw and not an entirely vector-only document.

The three saved curve files contain `roc_fpr`, `roc_tpr`, `pr_plot_recall`, and `pr_plot_precision`. PR AUC is HiLearn's native trapezoidal area including its artificial endpoints, not average precision. The script validates all three AUROC and PR AUC values before export.

## Provenance and validation

`provenance/original_reproduce_layout_v2.py` is the unmodified original assembly script, retained as a historical record. Its old absolute paths are not the portable entry point. `provenance/source_manifest.json` records hashes of bundled exact source assets, curves, and original deliverables. The supplied preview PDFs are the unchanged user-selected originals.

The portable default rebuild was compared with the original A4 PDF using the same renderer at 160 dpi: zero differing color components and maximum pixel difference 0. It was also rendered with Poppler. All 14 standalone panels were visually checked. See `provenance/portable_layout_validation.json` for the recorded result.

The original figure's embedded Arial text remains 5.25-8 pt at its specified final size. The original source-PDF hashes are:

```text
source_figure1_workflow_model.pdf
1e351185c1ab85cc4eaa937dec7abd9131002c1560a6166f80b61ccb45269676

source_figure2_tracks_correlations.pdf
f2e3c877ec46e5f83db358aa1a0c8ab796e3d0dc531edae0e8d09de071013c21
```
