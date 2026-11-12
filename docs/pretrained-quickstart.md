# Run pretrained GROPlant on real DNA

This example loads our trained **strand_seiplant_v1** models and performs actual
DNA-to-score inference. The model version retains repeats and uses **14 / 5 / 3
SRX** for Arabidopsis / wheat / maize. It does not use the separate three-library
wheat model or the N50 experiment.

## Install for CPU prediction

Use Python 3.10 or 3.11 in a separate environment. From the downloaded repository
root, run:

```bash
python -m pip install --upgrade pip setuptools wheel
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r requirements-inference.txt
python -m pip install --no-deps -e .
```

No GPU, BAM files, full reference genome or retraining is needed for this example.
Each species checkpoint is approximately 310 MB; download only the species you
need. Prediction uses CPU float32 with a small batch by default. The model has
approximately 77 million trainable parameters, so allow adequate available RAM.

## Download and check one species

```bash
python scripts/download_pretrained.py --species arabidopsis_thaliana
python scripts/check_pretrained_example.py --species arabidopsis_thaliana
```

The first command downloads the checkpoint and verifies its SHA-256 before
installing it under `results/strand_seiplant_v1/<species>/model/best.pt`.
The second command reads the bundled real DNA sequences, runs the network and
compares both predicted strands with the archived reference predictions. It
must exit successfully and report a passing comparison.

The other species identifiers are `triticum_aestivum` and `zea_mays`:

```bash
python scripts/download_pretrained.py --species all
python scripts/check_pretrained_example.py --species all
```

The check requires `abs(actual - expected) <= 0.0002 + 0.001 * abs(expected)`
for every strand score (`atol=2e-4`, `rtol=1e-3`). Original archived GPU
predictions and CPU float32 inference can differ slightly with hardware,
numerical kernels and batch size. The original expected answers are unchanged.
See [the validation record](../publication/metadata/pretrained_validation.json)
for observed differences. Outputs are saved in `results/pretrained_example/`
as one predictions TSV and one check JSON per species. Rerunning into an
existing output requires `--force`, or use a different `--output-dir`.
These selected examples are a software check, not a new performance estimate.

## Predict your own 1,024-bp sequences

```bash
python scripts/predict_fasta.py --species arabidopsis_thaliana \
  --fasta examples/real_dna/arabidopsis_thaliana.fa \
  --output predictions.tsv
```

On Windows PowerShell, enter the command on one line. Replace `--fasta` with
your own FASTA path to predict new sequences. Every record must have a unique
identifier and exactly 1,024 A/C/G/T bases. Lowercase is accepted; ambiguous
bases and incorrect lengths are rejected. The encoding retains the original
A/T/C/G channel order. FASTA sequences are used in the orientation supplied;
for genomic plus/minus interpretation supply the forward-reference sequence.

The output contains `sequence_id`, `prediction_plus` and `prediction_minus`.
Each score is between zero and one. Scores summarize learned consensus-derived
window targets; they are not CPM, calibrated probabilities or per-base signal.
The FASTA command does not reverse-complement or scan longer sequences.

For full-reference scanning, the downloaded checkpoints use the location
expected by `scripts/predict_genome.py`. Follow the reference preparation and
configuration instructions in [the model guide](publication-model.md).

## What is downloaded and how it is identified

The files are hosted as assets of the
[pretrained model release](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/strand-seiplant-v1-pretrained-20261112).
The small [download manifest](../publication/metadata/pretrained_models.json)
records each URL, file size, checksum and source checkpoint identity.

The release checkpoints retain all original model tensors and inference
metadata. The production configuration, which contained internal filesystem
paths, is omitted. Therefore the downloadable file checksum differs from the
original checkpoint checksum in `model_identity.json`; both are recorded, and
every model tensor is verified to be identical. The network architecture and
trained values are unchanged. Checkpoints are loaded with PyTorch's restricted
`weights_only=True` loader.

The actual DNA records, genomic coordinates, source identities, selection rule
and archived reference answers are under [examples/real_dna](../examples/real_dna/).
The full 193,508-window plotting and evaluation inputs remain under
`publication/data/`; those are separate from this small inference example.
