# GROPlant model

This folder contains the model implementation, parameters, pretrained-weight
download instructions and real DNA tests. The selected scientific model is
**`strand_seiplant_v1`**, repeats retained, with independent Arabidopsis / wheat /
maize models trained from **14 / 5 / 3 SRX**. Reorganizing the repository does
not change the trained weights or the archived prediction answers.

## Where to look

| Content | Location |
| --- | --- |
| Model architecture, alignment and label code | [src/plant_pausing_bowtie2/](src/plant_pausing_bowtie2/) |
| Selected model and label parameters | [configs/](configs/) |
| Training, preparation and prediction commands | [scripts/](scripts/) |
| Three weight URLs, file sizes and SHA-256 values | [metadata/pretrained_models.json](metadata/pretrained_models.json) |
| Model identity, assemblies and sample list | [metadata/](metadata/) |
| 72 real DNA examples, coordinates and expected predictions | [examples/real_dna/](examples/real_dna/) |
| Exclusion masks | [data/masks/](data/masks/) |
| Automated model and workflow checks | [tests/](tests/) |
| CPU inference guide | [docs/pretrained-quickstart.md](docs/pretrained-quickstart.md) |
| Architecture, labels and training guide | [docs/publication-model.md](docs/publication-model.md) |

Downloaded checkpoints and new predictions are written under `model/results/`.
The large checkpoint files are distributed through the
[fixed pretrained release](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/strand-seiplant-v1-pretrained-20261112),
rather than stored in Git. Reference genomes and original experimental inputs
for retraining must be obtained separately; configuration-relative paths point
to `model/references/` and `model/data/` by default.

## Run a real DNA test

Start with the [repository installation instructions](../README.md#start-from-a-fresh-download).
Keep the terminal at the **repository root**, then run:

```bash
python model/scripts/download_pretrained.py --species arabidopsis_thaliana
python model/scripts/check_pretrained_example.py --species arabidopsis_thaliana
```

Use `--species all` for all three species. Each successful species check reports
`passed: true`, 24 DNA sequences and 48 strand scores. Expected answers stay
unchanged; comparisons allow the documented floating-point tolerance.
These selected examples check software functionality, not independent model
performance. See the [CPU guide](docs/pretrained-quickstart.md) to predict your
own 1,024-bp FASTA sequences.

The full combined check is run from the repository root:

```bash
python check_installation.py --font "DejaVu Sans" --require-isolated-env
```

It also checks and redraws the results in [figure/](../figure/README.md).
Original code uses MIT; original weights and project-produced data use CC BY
4.0, with attribution **GROPlant contributors**. See the shared
[license scope](../docs/licensing.md) for third-party material and notices.
