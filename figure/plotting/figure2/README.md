# Figure 2: cross-species radar and genomic tracks

Reproduce the selected **2026-11-16 final A4 figure**, including its caption,
or export each of its six panels separately. All outputs are drawn from the
supplied numerical data as vector graphics; no image tracing is involved.

[Final A4 PDF](../../preview/figure2/GROseq_radar_tracks_with_legend_A4.pdf) ·
[Figure without caption](../../preview/figure2/GROseq_radar_tracks_artwork_183x210mm.pdf) ·
[Plotting data and definitions](../../data/figure2/README.md)

## Requirements

R 4.2.1 or later, the standard `grid` package, Cairo support, and installed
Arial are sufficient. No CRAN packages, trained models, GPU, BAM files, reference
genome or network access are required to redraw these results. Fonts are not
distributed. Another font may be substituted by the system if Arial is absent;
such output is not a typography-identical reproduction.

## Reproduce the complete figure

From the repository root:

```sh
Rscript --vanilla figure/plotting/figure2/plot_radar_tracks_composite.R
```

This creates `figure/results/figure2/figures/` with the A4 PDF and preview,
plus the 183 x 210 mm artwork PDF and 600-dpi PNG. It also writes numerical
validation tables under `figure/results/figure2/validation/`.
The original reference figures remain unchanged in `figure/preview/figure2/`.

Paths resolve relative to the script, so an absolute script path also works
from another working directory. Optional `--data-dir PATH` and `--output-dir PATH`
select alternate input and output folders. Repeating a command replaces only
its generated files in that output folder. The `--define-only` option loads
the data and drawing functions without creating figures.

## Export components or a single panel

```sh
# All six panels, two panel groups, and the shared legend (9 PDFs and PNGs)
Rscript --vanilla figure/plotting/figure2/plot_components.R --component all

# Only the three radar panels, or only the three track panels
Rscript --vanilla figure/plotting/figure2/plot_components.R --component radar
Rscript --vanilla figure/plotting/figure2/plot_components.R --component tracks

# One panel: replace a with b, c, d, e, or f
Rscript --vanilla figure/plotting/figure2/plot_components.R --component a
```

Component files are written to `figure/results/figure2/components/`.
`--component legend` exports the shared radar legend. Individual panels retain
their original letters, geometry, axis scales and data; they are translated to
their own pages without scaling. Single radar panels omit the shared legend,
which can be placed alongside them when assembling a new figure.

| Panel | Input species | Drawing function | Input data |
| --- | --- | --- | --- |
| a | Arabidopsis | `draw_radar(panels=1)` | The three metric tables, Arabidopsis rows |
| b | Maize | `draw_radar(panels=2)` | The three metric tables, maize rows |
| c | Wheat | `draw_radar(panels=3)` | The three metric tables, wheat rows |
| d | Arabidopsis | `draw_tracks(panels=1)` | Chr3:3,700,000-3,800,000, plus strand |
| e | Maize | `draw_tracks(panels=2)` | Chr6:119,800,000-119,900,000, minus strand |
| f | Wheat | `draw_tracks(panels=3)` | Chr3A:675,800,000-675,900,000, plus strand |

Both entry points use the functions in `plot_radar_tracks_composite.R`, so
the individual panels and complete composite share one implementation.
Detailed field mappings are in [provenance/panel_map.csv](provenance/panel_map.csv).
The editable caption is [figure_legend_submission.txt](figure_legend_submission.txt);
the palette is documented in [palette.json](palette.json), with the actual
colors defined by `SAME` and `CROSS` in the R script.

## Verify the inputs and regenerated tables

Python 3.10 or later, standard library only:

```sh
python figure/scripts/verify_radar_tracks_inputs.py
python figure/scripts/verify_radar_tracks_inputs.py --rendered-dir figure/results/figure2
```

The verifier checks all eight original data-file SHA-256 values, the nine
species/model combinations, pooled sample sizes, Jaccard from its stored counts,
and all nine regional correlations from the supplied track arrays. The second
command also checks regenerated numerical tables. Global Pearson and AUROC are
read from the frozen summaries; this is not a recalculation from raw reads or
the complete test-set prediction arrays.

The [reproduction record](provenance/reproduction_check.json) reports the actual
Windows / R 4.2.1 / Cairo / Arial check. The two composite PDFs rendered identically
to the selected originals at 144 dpi. All 11 generated PDFs were checked for page
geometry, embedded Arial, vector content and text bounds; the component renders
were visually inspected. This does not claim identical font rendering on every OS.

The exact original R script is retained in `provenance/` as a historical source.
Its former directory layout differs from the repository; use the two entry points
above for reproduction. The historical delivery manifest describes the original
delivery, while `figure/data/figure2/data_manifest.json` verifies current inputs.

After relocation into `figure2/`, both entry points were executed again with
their default input/output paths. The A4 and artwork PDFs matched the selected
originals pixel for pixel at 144 dpi; six individual panels, two panel groups
and one legend were exported. All 27 metric values and nine regional
correlations passed the numerical verifier. See the
[relocation validation](reference_validation/relocation_reproduction_20261118.json).

## 中文快速说明

本目录对应最终采用的雷达图与 tracks 组图。先运行整图命令，再按需运行
`plot_components.R --component a`（a-f 任一面板）即可单独导出。实际绘图数据在
`figure/data/figure2/`，原始最终 PDF 在 `figure/preview/figure2/`，
重画结果在 `figure/results/figure2/`。不需要重新训练或运行模型。

代码沿用仓库 MIT 许可；本项目生成的数据与图件沿用 CC BY 4.0，署名
GROPlant contributors，参见[许可说明](../../../docs/licensing.md)。
