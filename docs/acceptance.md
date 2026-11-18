# Installation and reproduction checks for GROPlant

`check_installation.py` checks pretrained CPU inference and bundled figure inputs. Reproducing the four selected final figures is optional (`--figures`). Neither mode trains a model or processes original FASTQ/BAM collections. Use a fresh clone/archive and record its source revision.

## Model environment and default check

Use a separate Python 3.10 or 3.11 environment:

```bash
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS or `.venv\Scripts\Activate.ps1` in Windows PowerShell. Do not use `--system-site-packages` for the isolated check. For Windows/Linux CPU PyTorch:

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r model/requirements-inference.txt -c reproducibility/requirements-acceptance-lock.txt
python -m pip install --no-deps -e .
python -m pip check
python check_installation.py --require-isolated-env
```

The unchanged lock records the historical Python 3.10/Linux CPU environment, not the four final-figure workflows. The installation commands target Windows/Linux; historical acceptance does not claim native macOS validation.

The default check requires no R, Cairo, Arial or PDF plotting packages. It stops at the first failed step:

1. Dependency consistency (`pip check`) and model unit tests.
2. Published file identities and consistency of all 193,508 test windows.
3. Selected final input/reference identities, Figure 2 radar/track checks, Figure S1's 48 regional summaries, and Figure S2's 60 nonzero decile estimates and saved bootstrap intervals.
4. Download or verify the three published checkpoints.
5. Predict all 72 real DNA sequences (144 strand scores) on CPU and compare with historical answers.

A first run downloads approximately 930 MB of checkpoints. Cached files are reused only after size/SHA-256 checks. Downloads use checked 4-MiB HTTP ranges with up to three attempts per range; complete-file responses to the first request are also supported. Failed downloads do not replace existing checkpoints. Proxy settings belong to the user's environment.

Score comparisons use `abs(actual-expected) <= 2e-4 + 1e-3*abs(expected)`. These examples check installation functionality, not representative model performance. Saved decile bootstrap draws are verified rather than resampled.

Require exit code zero and `"passed": true` in `reproducibility/acceptance/summary.json`. Its `final_figures` entry explicitly says `not run` without `--figures`. R is never searched for in default mode.

## Separate final-figure environment

| Figure | Content | Module |
| --- | --- | --- |
| Figure 1 | Model/workflow, tracks, correlations and ROC/PR | `figure/plotting/figure1/` |
| Figure 2 | Cross-species radar metrics and tracks | `figure/plotting/figure2/` |
| Figure S1 | Regional signal and correlations | `figure/plotting/figureS1/` |
| Figure S2 | CPM decile correlations, grouped strands | `figure/plotting/figureS2/` |

Use a separate Python 3.12 plotting environment: `python3.12 -m venv .venv-figures` on Linux, or `py -3.12 -m venv .venv-figures` on Windows. Activate it, then install existing module requirements together:

```bash
python -m pip install -r figure/plotting/figure1/requirements-source.txt -r figure/plotting/figure1/requirements-layout.txt -r figure/plotting/figureS1/requirements.txt -r figure/plotting/figureS2/requirements.txt
python -m pip check
```

Do not apply the historical acceptance lock here: `pypdf==7.0.1` and `reportlab==5.1.0` conflict with the final scripts' `pypdf>=6,<7` and `reportlab>=4,<5`. PyMuPDF and pdfplumber are also needed. These module requirements are not a newly certified cross-platform lock. Each reproduction records the actual plotting interpreter/package versions and separately runs its `pip check`.

Figure 2 requires R 4.2.1 or later with Cairo and standard R packages. In a shell supporting single quotes, check with `Rscript --vanilla -e 'stopifnot(capabilities("cairo"))'`. Supply `--rscript` when Rscript is not on PATH.

All four final figures require licensed Arial. Provide `arial.ttf`, `arialbd.ttf`, `ariali.ttf` and `arialbi.ttf` via `--font-dir` or `GROSEQ_FONT_DIR`. Arial must also be installed and discoverable by R/Cairo: a Python font directory does not install fonts for R. Font binaries are not distributed. Final A4 checks reject font substitution.

## Reproduce only the four final figures

From the plotting environment:

```bash
python figure/plotting/reproduce_all.py --rscript /absolute/path/to/Rscript --font-dir /absolute/path/to/Arial
```

This entry performs no model download or inference. It can also be launched from another environment with `--python /absolute/path/to/plotting/python`. Inputs resolve from the script location, so launch from outside the repository is supported. `--dry-run` prints commands without execution or file creation.

Figure 1 first redraws its numerical source tracks/scatters, then passes the generated source PDF to the assembler. Original vector assets remain authoritative for workflow/model panels; ROC/PR uses saved vertices. The runner then produces Figure 2 and its components, Figure S1, and Figure S2 and its components. Figure 2's rendered values are checked against input summaries. Builders perform their numerical/layout checks. The final check requires four nonempty, single-page A4 PDFs with extractable Arial text. Structural checks do not replace visual review or claim pixel equivalence.

Raw-input HiLearn ROC/PR recalculation remains available under `figure/plotting/figure1/roc_pr/plot_publication_roc_pr.py`, with additional HiLearn dependencies described there. It is not a fifth final-figure workflow and is not automatically run.

## Model check plus all four figures

From the model environment, select the plotting interpreter:

```bash
python check_installation.py --require-isolated-env --figures --figure-python /absolute/path/to/plotting/python --rscript /absolute/path/to/Rscript --font-dir /absolute/path/to/Arial
```

On Windows use quoted absolute `.exe` paths. The legacy `--font` option remains accepted by older model-only commands; `--figures` requires Arial. Historical DejaVu Sans results establish numerical reproduction, not these final figures' typography.

The combined summary embeds the figure report only after all steps succeed. Model-only success is not labelled complete figure reproduction. `--require-isolated-env` checks the model interpreter; the plotting interpreter is separately identified in its report.

## Output locations and historical evidence

| Output | Location |
| --- | --- |
| Checkpoints | `model/results/strand_seiplant_v1/<species>/model/best.pt` |
| Predictions and per-species checks | `model/results/pretrained_example/` |
| Figure 1 and regenerated source panels | `figure/results/figure1/` |
| Figure 2 and components | `figure/results/figure2/` |
| Figure S1 and components | `figure/results/figureS1/` |
| Figure S2 and components | `figure/results/figureS2/` |
| Standalone figure logs/summary | `figure/results/reproduction_validation/` |
| Model acceptance logs/summary | `reproducibility/acceptance/` |
| Figure logs for combined check | `reproducibility/acceptance/final_figures/` |

The figure runner accepts `--output-dir` and `--report-dir`. The model checker accepts `--output-dir` for its acceptance logs.

The [v0.4.0 release](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/groplant-v0.4.0-20261112) and [v0.3.3 release](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/groplant-v0.3.3-20261112) preserve historical installation/download/inference evidence for their exact revisions. They do not certify this later four-figure runner or its separate plotting environment. Current-run logs establish which checks actually ran on the current source.
