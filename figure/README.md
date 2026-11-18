# GROPlant figures and results

## Selected final figures

The manuscript figure number is also the directory name under `plotting/`,
`data/` and `preview/`. Start with the row for the figure you want to reproduce.

| Figure | Content | Drawing scripts and commands | Plot data | Original PDF |
|---|---|---|---|---|
| **Figure 1** | Workflow, model, genomic tracks, scatter and ROC/PR; a–n | [figure1/](plotting/figure1/README.md) | [figure1/](data/figure1/) | [A4 PDF](preview/figure1/GROseq_combined_v2_A4.pdf) |
| **Figure 2** | Cross-species radar and genomic tracks; a–f | [figure2/](plotting/figure2/README.md) | [figure2/](data/figure2/) | [A4 PDF](preview/figure2/GROseq_radar_tracks_with_legend_A4.pdf) |
| **Figure S1** | Regional correlation and mean experimental signal; a–f | [figureS1/](plotting/figureS1/README.md) | [summary](data/figureS1/) + [shared windows](data/regional/) | [A4 PDF](preview/figureS1/region_signal_six_panels_A4.pdf) |
| **Figure S2** | CPM deciles with separately estimated strands; a–c | [figureS2/](plotting/figureS2/README.md) | [figureS2/](data/figureS2/) | [A4 PDF](preview/figureS2/reference/GROseq_CPM_deciles_grouped_strands_A4.pdf) |

[Complete figure previews](preview/README.md) ·
[Entry filenames, run commands and panel export](plotting/README.md#selected-final-figures-script-index)

The original PDF filenames are retained for provenance; their folder names
identify the manuscript numbers. Each entry exports the full figure and separate
panels. Figure 1 a/b use the exact author-provided vector PDF; the complete
original program that created those schematic primitives was not recovered.

## Folder guide

| Location | Contents |
|---|---|
| `plotting/figure1, figure2, figureS1, figureS2/` | The four final drawing modules, settings and panel mappings |
| `data/figure1, figure2, figureS1, figureS2/` | Corresponding frozen plotting inputs |
| `data/regional/` | Shared full-window data for S1/S2 and consistency checks; not an extra figure |
| `preview/figure1, figure2, figureS1, figureS2/` | Original PDFs, complete previews and reference records |
| `scripts/` | Input-verification tools only |
| `docs/` | Quantity definitions and statistical conventions |
| `reproducibility/` | Recorded software environment |
| `results/` (local only) | Regenerated outputs, excluded from Git |

Figure 1 also retains [ROC/PR reconstruction](plotting/figure1/roc_pr/README.md)
from [complete truth labels and scores](data/figure1/roc_inputs/). This is a
supporting numerical check within Figure 1, not a fifth figure.

These results use **`strand_seiplant_v1`**, repeats retained, Arabidopsis /
wheat / maize = **14 / 5 / 3 SRX**, covering **193,508** held-out test windows.
Plotting does not need model weights, FASTQ/BAM files, a genome FASTA or a GPU.

## Reproduce

Use the [plotting guide](plotting/README.md) for the exact commands and a separate
Python 3.12 plotting environment. Figure 2 also uses standard R with Cairo.
Install the module requirements and locally available Arial fonts as documented;
font binaries are not redistributed.

To draw all four final figures after setup, run from the repository root:

```sh
python figure/plotting/reproduce_all.py --font-dir "C:/Windows/Fonts"
```

Use `--rscript` if Rscript is not on PATH. Individual entries remain available
in each module README. Generated files go to `figure/results/`, leaving the
versioned inputs and original previews unchanged.

## Check the data

```sh
python figure/scripts/verify_publication_inputs.py
python figure/scripts/verify_selected_final_inputs.py
python figure/scripts/verify_radar_tracks_inputs.py
python figure/scripts/verify_decile_inputs.py
```

The checks verify file identities, full-window consistency, the radar/track
inputs and the decile estimates/intervals. Plotting S1 also independently checks
its regional summaries against every window. Frozen bootstrap draws are supplied;
redrawing S2 does not resample them.

The [installation guide](../docs/acceptance.md) distinguishes model validation
from optional figure reproduction. The [figure methods](docs/publication-figures.md)
define the plotted quantities.

Code uses MIT. Project-produced data and figures use CC BY 4.0 with attribution
**GROPlant contributors**; see [license scope](../docs/licensing.md).
