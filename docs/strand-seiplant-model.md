# Species-specific DNA-to-GRO-seq strand model

This workflow trains one model per species to predict the strand-resolved GRO-seq labels
`y_plus` and `y_minus` directly from a 1,024-bp genomic DNA sequence. It adapts the
SeiPlant architecture and genome-scanning strategy described by Lv *et al.* (2026),
*Cross-species prediction of histone modifications in plants via deep learning*
([DOI: 10.1186/s13059-025-03929-4](https://doi.org/10.1186/s13059-025-03929-4)).
The reference implementation consulted during development was
[`compbioNJU/SeiPlant`](https://github.com/compbioNJU/SeiPlant), commit
`1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853`.

The model predicts the previously constructed GRO-seq interval-coverage labels. These
targets measure reproducible nascent-transcription coverage and are not single-base pause
sites, read counts, expression estimates, or probabilities.

## What is retained from SeiPlant

- 1,024-bp forward-reference DNA input encoded in four channels ordered A, T, C, G;
- three hierarchical convolution stages with 480, 640, and 960 channels;
- kernel size 9 in the local and residual blocks, max pooling by 4 after stages 1 and 2;
- dropout 0.2 in stages 2 and 3;
- five residual dilated convolutions with dilation 2, 4, 8, 16, and 25, kernel size 5,
  960 channels, and dropout 0.1;
- compression to 16 cubic B-spline basis functions;
- sigmoid-bounded output and MSE regression loss;
- Adam optimization with learning rate `1e-5`;
- chromosome-separated train, validation, and test sets;
- 1,024-bp genome inference windows at 128-bp steps, with only positions 448–576
  assigned the window prediction.

The article calls the B-spline transformation learnable, whereas the cited public
SeiPlant implementation fixes the spline basis and trains the subsequent layers. This
repository follows the public implementation and records that choice explicitly.

## Project-specific adaptations

1. The output has two continuous targets: plus-strand and minus-strand GRO-seq signal.
   Both are learned jointly from a shared sequence backbone, but have distinct output
   weights. They are never summed or collapsed.
2. A 256-unit ReLU head is used before the two output units. The original code used the
   number of chromatin targets as the hidden width; using that literal rule here would
   leave only two hidden units.
3. The released y tables contain only windows where at least one strand is positive. A
   deterministic set of zero-signal genomic windows is therefore sampled separately in
   each chromosome split. The default background:positive ratio is 1:1 and can be changed
   in the configuration. This is necessary to calibrate genome-wide background output.
4. Raw predictions are retained. The paper's optional `<0.01` threshold and subsequent
   min-max display normalization are not applied to quantitative output because they
   would change the model scale. A display-only transformation can be added downstream.
5. The example configuration follows the article's strict A/T/C/G-only filter. The
   `max_ambiguous_fraction` option can be increased for a documented sensitivity analysis;
   retained ambiguous positions are encoded as all-zero columns.

## Chromosome splits

| Species | Training | Validation | Test | Rationale |
|---|---|---|---|---|
| *Arabidopsis thaliana* | 1, 2, 5 | 4 | 3 | Same validation/test chromosomes as the article |
| *Zea mays* | 1–5, 7, 9, 10 | 8 | 6 | Same validation/test chromosomes as the article |
| *Triticum aestivum* | homoeologous groups 1, 2, 5, 6, 7 | 4A/4B/4D | 3A/3B/3D | Article-analog extension that keeps all A/B/D homoeologs in one split |

Organellar sequences, unplaced wheat sequence `Un`, and maize scaffolds are excluded.
Keeping wheat homoeologs together reduces information leakage caused by highly similar
A-, B-, and D-subgenome sequence.

## Install and inspect

The GPU stack is separate from the alignment environment:

```bash
conda env create -f environment-model.yml
conda activate plant-pausing-model
python -m pip install --no-deps -e .
plant-pausing-model plan configs/example_strand_model.json
```

Relative paths are resolved from the JSON file. Copy the example into the analysis
project's `configs/` directory if the reference genomes and results live outside the Git
checkout.

## Prepare manifests

```bash
python scripts/prepare_model_data.py configs/example_strand_model.json
```

For each species, preparation performs the following operations independently in the
train, validation, and test chromosome sets:

1. read the released positive y table and retain only configured nuclear chromosomes;
2. remove windows overlapping the species mask or exceeding the ambiguous-base limit;
3. sample background windows from the same 1,024/512-bp genomic grid without replacing
   positive coordinates;
4. order retained windows by chromosome and coordinate and encode each base once into a
   memory-mapped `uint8` sequence cache;
5. write `train.tsv.gz`, `validation.tsv.gz`, and `test.tsv.gz` plus input hashes and counts.

The background sample never crosses chromosome splits. Source is recorded as `positive`
or `background`; background rows have `y_plus=y_minus=0`.

## Train one model per species

```bash
python scripts/train_strand_model.py configs/example_strand_model.json \
  --species arabidopsis_thaliana --device cuda:0
python scripts/train_strand_model.py configs/example_strand_model.json \
  --species triticum_aestivum --device cuda:1
python scripts/train_strand_model.py configs/example_strand_model.json \
  --species zea_mays --device cuda:2
```

The three commands can run concurrently on different GPUs. Training uses mixed precision
when CUDA is available, a batch size of 256, up to 30 epochs, and early stopping after 5
epochs without validation-MSE improvement. The best checkpoint is selected only by the
validation chromosomes. The held-out test chromosomes are evaluated once after model
selection. Reported metrics are MSE, MAE, Pearson correlation, and Spearman correlation,
separately for plus and minus strands.

On the production RTX 3080 Ti, a real forward/backward smoke test at batch size 256 used
approximately 3.2 GiB peak allocated GPU memory. This leaves headroom for CUDA workspace
and makes three concurrent species-specific jobs practical on separate GPUs.

Each species produces:

```text
results/strand_seiplant_v1/<species>/
  manifests/
    train.tsv.gz
    validation.tsv.gz
    test.tsv.gz
    train.sequences.uint8
    validation.sequences.uint8
    test.sequences.uint8
    summary.json
  model/
    best.pt
    history.tsv
    validation_predictions.npz
    test_predictions.npz
    metrics.json
```

Each cache stores one byte per base (`A=0`, `T=1`, `C=2`, `G=3`, other=4) in exactly
the same row order as its manifest. It is memory mapped by data-loader workers, avoiding
hundreds of thousands of random FASTA seeks in every epoch. Training batches are still
shuffled deterministically by the PyTorch generator seeded with 42.

## Genome-wide prediction

```bash
python scripts/predict_genome.py configs/example_strand_model.json \
  --species arabidopsis_thaliana --device cuda:0
```

The command scans the configured nuclear chromosomes in 1,024-bp windows with a 128-bp
step. A scalar strand score predicted from each sequence is assigned to its central 128 bp.
It writes three raw BedGraph files:

- `*.plus.bedGraph`: positive plus-strand predictions;
- `*.minus.bedGraph`: positive minus-strand predictions;
- `*.minus_signed.bedGraph`: the same minus values multiplied by -1 for mirrored plots.

The minus-signed file is a visualization representation only. Statistical evaluation uses
the original non-negative `minus` values.

For a rapid browser-style figure before whole-genome inference, restrict prediction to a
bounded region; commas in coordinates are optional:

```bash
python scripts/predict_genome.py configs/example_strand_model.json \
  --species arabidopsis_thaliana --device cuda:0 --force \
  --region 3:9141194-9207074
```

The region uses the same 1,024/128-bp inference and central-bin rule as a complete scan.
It writes the standard track names, so a later chromosome or whole-genome run must use
`--force` to replace this bounded result.

## Plot observed and predicted tracks

```bash
python scripts/plot_strand_tracks.py configs/example_strand_model.json \
  --species arabidopsis_thaliana --chrom 3 \
  --start 9000000 --end 9100000 \
  --output results/strand_seiplant_v1/figures/arabidopsis_chr3
```

The output PNG and vector PDF contain observed and predicted panels. Plus-strand signal
is plotted above zero and minus-strand signal below zero. The user's example corresponds
to the track view in Figure 1c of Lv *et al.*; Figure 3c in that article is a different
family-level comparison panel.

## Reproducibility and interpretation

- `--force` is required to replace manifests, checkpoints, or complete prediction output.
- The manifest summary hashes the source y file, FASTA index, and generated manifests.
- The model metadata records the PyTorch version, seed, best epoch, checkpoint hash, and
  strand-specific validation/test metrics.
- Reference FASTA and y files are intentionally excluded from Git; their release paths and
  hashes belong in the project result metadata.
- A high predicted score means sequence context resembles reproducible strand-specific
  nascent-transcription intervals in the training data. It does not by itself establish a
  pause site or enhancer.
