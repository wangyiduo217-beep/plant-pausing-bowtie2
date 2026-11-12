# Small real-DNA prediction examples

These fixtures let a new user load this project's pretrained GRO-seq weights,
predict actual plant DNA, and compare the two output scores with the frozen
publication results. See the repository's [quick start](../../../README.md) for
installation, weight download, prediction, and verification commands.

## Included files

| Species | FASTA | Reference predictions | Test chromosomes | Sequences |
| --- | --- | --- | --- | ---: |
| *Arabidopsis thaliana* | `arabidopsis_thaliana.fa` | `arabidopsis_thaliana.expected.tsv` | 3 | 24 |
| *Triticum aestivum* | `triticum_aestivum.fa` | `triticum_aestivum.expected.tsv` | 3A, 3B, 3D; 8 per chromosome | 24 |
| *Zea mays* | `zea_mays.fa` | `zea_mays.expected.tsv` | 6 | 24 |

Every FASTA record is one **1,024-bp, A/T/C/G-only, forward-reference sequence**
from the original held-out test set. No full reference genome, BAM, or FASTA
index is needed to run this small FASTA prediction example. The three pretrained
models are separate downloads. Use each species' FASTA with its matching model.

The experiment is **`strand_seiplant_v1`, repeats retained, 14/5/3 SRX**. In
particular, these wheat sequences and answers belong to the original five-SRX
model, not to the separate three-SRX retraining or N50 experiment.

## What the reference table contains

- `sequence_id`: the first token in the matching FASTA header; unique within
  this package, formatted as `<species_slug>_w<window_index>`.
- `prediction_plus`, `prediction_minus`: the original, unmodified **raw model
  scores** in the committed publication plotting table. These are not CPM,
  calibrated probabilities, or base-resolution measurements. They are the
  values to compare against freshly generated model predictions.
- `window_index`: the original one-based row number in the species' test
  manifest, excluding its header.
- `chrom`, `start0`, `end0`: reference chromosome and 0-based, half-open genomic
  coordinates. `end0 - start0` is always 1,024.
- `sequence_sha256`: SHA-256 of the uppercase 1,024-base sequence alone,
  excluding the FASTA header and line breaks.
- `selection_rank0`, `selection_pool_size`: the selected rank and available
  window count within that test chromosome, under the rule below.
- `y_plus`, `y_minus`, `label_source`: the archived consensus-derived training
  targets and the original `positive`/`background` assignment. These provide
  context and are **not** the reference predictions. They must not be confused
  with experimental CPM or used as exact expected network output.

Floating-point implementations can yield small numerical differences. Use the
documented checker tolerance, rather than demanding byte-identical score text.

## Deterministic selection and provenance

Within each held-out chromosome, windows are ordered by the larger of their
two frozen prediction scores, breaking ties by integer `window_index`. For a
chromosome containing `n` windows, select ranks
`floor(i * (n - 1) / (k - 1))`, for `i = 0, ..., k - 1`. Here `k = 24` for
Arabidopsis and maize and `k = 8` for each of the three wheat test chromosomes.
This deliberately spans low and high scores, including extremes. No random
seed or model retraining is involved.

The source prediction tables are the committed
[`figure/data/regional/`](../../../figure/data/regional/) tables. DNA
was decoded from the corresponding original production test sequence caches
(`A=0, T=1, C=2, G=3`) after verifying each complete cache and test-manifest
SHA-256 against the published
[`model_identity.json`](../../metadata/model_identity.json).
Every selected sequence was then independently extracted by coordinate from
the matching reference FASTA and required to match exactly. The reference
FASTA-index hashes were also verified. No production input was changed.

[`manifest.json`](manifest.json) records all selected row indices, source-table
and fixture hashes, verified manifest/cache/index hashes, reference assembly
accessions and download URLs, and the required original checkpoint hashes.
The assemblies are TAIR10, IWGSC, and Zm-B73-REFERENCE-NAM-5.0, respectively.

**These are functionality examples, not a representative performance benchmark.**
They are intentionally selected using existing model scores. Do not estimate
model accuracy from these 24 windows or describe them as independent biological
replicates. The original full held-out test sets remain available as publication
plotting/evaluation inputs.
