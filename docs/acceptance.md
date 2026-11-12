# Installation and reproduction checks for GROPlant

This check covers the workflows distributed with the repository: pretrained
CPU inference, all six R figure families, HiLearn ROC/PR, and the identities and
scientific values of the bundled inputs. It does not retrain the models or
process the original FASTQ/BAM collections. All commands below run from the
repository root and use the `model/` and `figure/` layout introduced in 0.4.0.

## 1. Obtain a new copy

Use **Code → Download ZIP** on GitHub and extract the archive into a new folder,
or clone the repository into an empty destination. Record the downloaded
commit or release. The release validation report identifies the exact tested
commit and tree; a later `main` can contain later changes.

## 2. Install an independent Python environment

Use Python 3.10 or 3.11 in the new repository folder:

```bash
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS, or
`.venv\Scripts\Activate.ps1` in Windows PowerShell. Do not add
`--system-site-packages` when recreating the isolated check.

For Windows/Linux CPU PyTorch:

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r model/requirements-inference.txt -r figure/requirements.txt -c reproducibility/requirements-acceptance-lock.txt
python -m pip install --no-deps -e .
python -m pip check
```

The lock file records the tested Python 3.10/Linux CPU dependency versions.
`pip check` should report `No broken requirements found.` The installation
guide above is intended for Windows/Linux; a native macOS PyTorch install
uses the appropriate official wheel rather than the Linux/Windows CPU index.
The published acceptance record does not claim a macOS validation.

## 3. Provide R and a font

Use R 4.2.1 or later with Cairo enabled. Only standard R packages are needed.
`Rscript --vanilla -e 'stopifnot(capabilities("cairo"))'` checks Cairo in a
shell that supports single quoting. In R itself, run
`stopifnot(capabilities("cairo"))`.

The independent Linux validation installed a separate R 4.4.3 environment.
One way to install a separate R runtime is:

```bash
conda create -n groplant-r -c conda-forge r-base=4.4
```

Keep the Python venv active for the Python commands and pass the R environment's
`bin/Rscript` to `--rscript`. Adding that R `bin` directory to PATH also makes
its font-discovery tools available to Python plotting.

`figure/reproducibility/R-linux-64-explicit.txt` records its exact package URLs for
Linux; it is not a Windows/macOS environment file.

Arial is the requested publication font. The independent Linux acceptance
uses installed **DejaVu Sans**, explicitly passed to both Python and R; these
checks establish numerical reproduction and working output generation, not
pixel-identical Arial typography. The original Windows R previews were also
rendered and inspected. Font binaries are not bundled in the repository.

## 4. Run the complete check

```bash
python check_installation.py --font "DejaVu Sans" --require-isolated-env
```

If Rscript is not on PATH:

```bash
python check_installation.py --rscript "/absolute/path/to/Rscript" --font "DejaVu Sans" --require-isolated-env
```

For a Windows installation with Arial, use `--font "Arial"` and the actual
installed `Rscript.exe` path. A new run downloads three checkpoints totaling
about 930 MB. Existing checkpoints are reused only after size and SHA-256
checks. The output reports whether files were newly downloaded or verified.
Downloads use 4-MiB HTTP ranges, with at most three attempts per range. Each
range must match its requested position and size, and the completed file must
match the published SHA-256. A server that ignores the first range request can
return a complete file instead. Failed attempts never replace an existing
checkpoint; temporary partial files are removed when the command fails.
HTTP proxy settings, where needed, belong to the user's network configuration;
no credentials or proxy endpoints are embedded in the project.

The check runs these steps and stops at the first failure:

1. Installed dependency consistency (`pip check`).
2. Software unit tests.
3. File identities and all 193,508 formal test-window inputs.
4. All 60 nonzero decile estimates, saved bootstrap intervals and 288 frequency bins.
5. Download and verify the three pretrained checkpoint files.
6. Actually predict all 72 real DNA windows (144 strand scores) on CPU.
7. Check R/Cairo and generate the six R figure families.
8. Generate ROC/PR using the unchanged HiLearn 0.2.2 implementation.

Successful completion prints:

```text
PASS: all acceptance checks completed; 72 real DNA sequences / 144 strand scores checked.
```

The exit code must be zero. `reproducibility/acceptance/summary.json` must contain
`"passed": true`, and every recorded step must pass. Logs are saved beside
that summary. The isolated-environment flag is recorded in the summary.

## Output locations

| Output | Location |
| --- | --- |
| Checkpoints | `model/results/strand_seiplant_v1/<species>/model/best.pt` |
| Real DNA predictions and per-species comparison JSON | `model/results/pretrained_example/` |
| Six R figure families | `figure/results/correlations/`, `tracks/`, `regional_bars/`, `dual_axis_violins/`, `deciles/`, `frequency/` |
| ROC/PR figures and metrics | `figure/results/roc_pr/` |
| Complete acceptance summary and step logs | `reproducibility/acceptance/` |

The R output contains A4 proofs and 183-mm-wide vector PDFs plus 600-dpi PNGs.
The default decile layout groups both strands in three species panels; the
frequency module also creates a logarithmic-y detail view. The optional
six-panel decile layout is described in `figure/plotting/README.md` and
was checked separately.

CPU predictions are compared against unchanged historical answers using
`abs(actual-expected) <= 2e-4 + 1e-3*abs(expected)` for every strand score.
The examples check functionality and are not representative performance
estimates. The R decile plot reads frozen estimates and intervals; 1,000 saved
bootstrap replicates are supplied for verification, and no resampling is needed
to redraw the plot. The frequency plot
independently recomputes its bins from all raw window values.

## Verification record

The [GROPlant v0.4.0 release](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/groplant-v0.4.0-20261112)
provides the regression check for this reorganized layout, identifying its
exact source commit and output paths. This check reuses the independent
dependency environment and verified checkpoint cache established for
[v0.3.3](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/groplant-v0.3.3-20261112).
The earlier release preserves the evidence of fresh dependency installation
and public weight downloads. Reports distinguish new source downloads from
reused dependencies or verified weights, and record the actual font and
platform limits. The steps above also support a new user's fresh installation.
