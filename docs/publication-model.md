# Publication model: strand_seiplant_v1

This page specifies the model and input cohort underlying the current publication
figures. It freezes the **repeat-retained, 14 / 5 / 3 SRX version** for
*Arabidopsis thaliana*, *Triticum aestivum* and *Zea mays*, respectively. All five
wheat libraries are retained. The N50 experiment and the separate three-library
wheat experiment are outside this release.

The executable model and label implementations are unchanged:

- [`strand_model.py`](../src/plant_pausing_bowtie2/strand_model.py): preparation,
  architecture, training, evaluation and genome scanning;
- [`y_labels.py`](../src/plant_pausing_bowtie2/y_labels.py): per-library GRO-seq
  intervals, independent-library consensus and training targets;
- [`publication_v1_model.json`](../configs/publication_v1_model.json) and
  [`publication_v1_labels.json`](../configs/publication_v1_labels.json): connected,
  relative-path configurations for this version;
- [`model_identity.json`](../publication/metadata/model_identity.json): code,
  reference, label, manifest, sequence-cache and selected-checkpoint identities;
- [`samples.tsv`](../publication/metadata/samples.tsv): all 22 SRX and 34 SRR
  accessions, library layouts, grouping and strand rules.

## Inputs and reference assemblies

Each observation is one 1,024-bp sequence on the forward reference genome. Input
channels are ordered **A, T, C, G**, producing a `4 x 1024` one-hot array.
The publication configuration requires A/T/C/G-only windows
(`max_ambiguous_fraction = 0`). Upper/lower case is normalized. Genome repeat
annotations are **not** used to remove repeat-containing windows or to replace
repeat bases with N.

References originate from Ensembl Plants release 63:

| Species | Assembly | Assembly accession | Source FASTA |
|---|---|---|---|
| *A. thaliana* | TAIR10 | GCA_000001735.1 | [TAIR10 toplevel FASTA](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) |
| *T. aestivum* | IWGSC | GCA_900519105.1 | [IWGSC toplevel FASTA](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz) |
| *Z. mays* | Zm-B73-REFERENCE-NAM-5.0 | GCA_902167145.1 | [NAM-5.0 toplevel FASTA](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz) |

The wheat reference is the IWGSC assembly identified above, **not RefSeq v2.1**.
Use the exact assembly and sequence names in the configuration; do not substitute
a newer reference without treating that as a different experiment. Downloaded
FASTA-gzip checksums and the recorded FASTA-index checksums appear separately in
the identity file. They identify different files. After decompression, place each
reference at `references/<species_slug>/genome.fa` and build its index with
`samtools faidx`.

### rRNA/organellar exclusions are retained

Repeat retention does not remove the existing rRNA/organellar or ambiguous-base
restrictions. The three small archived coordinate masks are supplied under
[`publication/data/masks/`](../publication/data/masks/), and both configurations
point to those exact files:

- Arabidopsis: 10 BED records on chromosomes 2/3 and Mt/Pt, covering known
  rRNA-enriched regions and organellar sequence;
- wheat: 1,946 BED records representing the known rRNA-enriched coordinate mask;
- maize: an **intentionally empty, zero-byte BED file** in this archived version.
  Do not replace it with a repeat-annotation mask. The maize BAMs had already
  undergone the recorded sequence-level rRNA filtering.

Mask checksums, record counts and reference-sequence names are recorded in
`model_identity.json`. Coordinates use 0-based, half-open BED intervals.

## Constructing the two training targets

The model learns reproducible interval coverage, **not CPM or a calibrated class
probability**. The two output targets are `y_plus` and `y_minus` for the forward
reference plus and minus strands. Neither is an inferred gene-relative sense
label, and the two are not summed.

Start with the archived processing definition: primary MAPQ >= 20 alignments,
excluding SAM flags 4, 256, 512 and 2048 (`-F 2820`). Coordinate duplicates are
not removed. Single-end libraries contribute their eligible reads; paired-end
libraries contribute R1 only (`-f 64`). Multiple SRRs in one SRX form one
independent library before consensus calculation. GSE181488 has RNA strand
opposite to alignment strand and uses HOMER `-rev`; the other retained libraries
use the same strand. These rules are explicit per library in the configuration.

Per SRX, the label workflow calls HOMER with:

```bash
findPeaks <tag_directory> -style groseq -tssFold 4 -bodyFold 3 \
  -minBodySize 500 -maxBodySize 100000 -pseudoCount 1 -o auto
```

The workflow adds `-rev` only where configured. It then builds independent
species/strand consensus with `bedtools multiinter` and retains atomic intervals
supported by at least two independent SRX libraries. If `q_i = support_i / N`,
the confidence weight is:

```text
c_i = 0.1 + 0.9 * (q_i - q_min) / (q_max - q_min)
```

Here `q_min` and `q_max` are observed among retained intervals for the same
species and strand. If all observed support values are equal, `c_i = 1`.
Keep the implemented operation order, including division by N before scaling;
algebraically simplifying it can alter rare floating-point rounding boundaries.

For each 1,024-bp window on the fixed 512-bp grid and each strand:

```text
raw_score = sum(overlap_bp(window, consensus_interval_i) / 1024 * c_i)
y = min(floor(raw_score * 10 + 0.5) / 10, 1)
```

Values below 0.1 become zero. The positive-window table retains rows where at
least one rounded strand target is nonzero. Its `binary_plus/binary_minus`
columns instead indicate **raw overlap before rounding**; they need not equal
`y >= 0.1` on every retained row. Model preparation separately samples 1:1
zero-labelled background windows on the same genomic grid within each split,
excluding occupied positive coordinates, coordinate masks and ambiguous DNA.

The label settings and actual SRX/SRR inputs are fully specified in the
publication label configuration. Its output directory matches the `positive_y`
paths in the publication model configuration. Detailed algorithm documentation
is also available in [the label workflow](groseq-y-label-generation.md).

## Architecture and training

Dimensions below exclude the batch dimension and use channels x positions.

| Step | Output | Parameters and operation |
|---|---|---|
| DNA encoding | `4 x 1024` | A/T/C/G one-hot channels |
| Local residual stage | `480 x 1024` | Two linear Conv1d layers, then two Conv1d + ReLU layers; `k=9`, padding 4; add the linear-stage output |
| Intermediate stage 1 | `640 x 256` | MaxPool kernel/stride 4, dropout 0.2; the same two-linear/two-nonlinear residual pattern, `k=9` |
| Intermediate stage 2 | `960 x 64` | A second pool-by-4 and residual stage; dropout 0.2, `k=9` |
| Dilated residual blocks | `960 x 64` | Five sequential blocks; `k=5`, dilations 2/4/8/16/25, padding `2*d`, dropout 0.1; each adds `ReLU(conv(dropout(x)))` to x |
| B-spline aggregation | `960 x 16` | Fixed cubic, degree-3 basis with 16 basis functions; multiply the 64-position feature map by a `64 x 16` basis matrix |
| Shared regression head | `15360 -> 256 -> 2` | Flatten; dropout 0.5; Linear + ReLU; Linear + elementwise Sigmoid |

The B-spline basis is fixed, not learned. Subsequent trainable layers learn to
use the aggregated features. Both strands share the backbone/head and have
separate final output weights. Training minimizes mean squared error across the
two continuous outputs using Adam, learning rate `1e-5`, weight decay 0, batch
size 256 and seed 42. The schedule permits up to 30 epochs, with early stopping
after five consecutive epochs without validation-MSE improvement. CUDA automatic
mixed precision is enabled by the configuration.

### Chromosome splits and sample counts

| Species | Training chromosomes | Validation | Held-out test |
|---|---|---|---|
| Arabidopsis | 1, 2, 5 | 4 | 3 |
| Wheat | 1A/B/D, 2A/B/D, 5A/B/D, 6A/B/D, 7A/B/D | 4A, 4B, 4D | 3A, 3B, 3D |
| Maize | 1-5, 7, 9, 10 | 8 | 6 |

Wheat homoeologous chromosome groups stay together. Unplaced and organellar
sequences are outside these model splits. Archived totals include the sampled
background windows:

| Species | Train windows | Validation windows | Test windows | Independent SRX |
|---|---:|---:|---:|---:|
| Arabidopsis | 148,862 | 34,534 | 45,066 | 14 |
| Wheat | 454,892 | 86,412 | 98,656 | 5 |
| Maize | 500,228 | 52,736 | 49,786 | 3 |

Window counts are not biological replicate counts. Checkpoint selection uses
validation chromosomes; test chromosomes are evaluated after selection. The
archived best epochs are 12, 10 and 8, respectively. The production runs include
resume events: model weights were restored while Adam moments were reinitialized,
as recorded in the identity metadata. A fresh uninterrupted run is not claimed
to recreate the archived weights bit for bit.

## Reproduction commands

The core model uses Python, PyTorch, NumPy and SciPy. The supplied model environment
also includes Matplotlib. The archived production framework was
`PyTorch 2.6.0+cu124`; the environment specification gives dependency bounds,
not a fully locked copy of every production package. Label execution requires
Linux, HOMER, SAMtools and BEDTools. The processing environment and tool-version
records are provided separately. No training or large-data pipeline starts when
using a `plan` command.

```bash
conda env create -f environment-model.yml
conda activate plant-pausing-model
python -m pip install --no-deps -e .

# Inspect all paths and scientific settings without reading large inputs.
plant-groseq-y plan configs/publication_v1_labels.json
plant-pausing-model plan configs/publication_v1_model.json

# Use the separate environment containing the biological command-line tools.
conda env create -f environment.yml
conda activate plant-pausing-bowtie2
python -m pip install --no-deps -e .

# With the required reference indexes, masks and processed BAMs available:
python scripts/call_groseq_intervals.py configs/publication_v1_labels.json
python scripts/build_groseq_y.py configs/publication_v1_labels.json

# Switch back to the model environment for sequence preparation and training.
conda activate plant-pausing-model

# Prepare ALL THREE species together, in the supplied configuration order.
python scripts/prepare_model_data.py configs/publication_v1_model.json

# Species-specific training may then be launched separately.
python scripts/train_strand_model.py configs/publication_v1_model.json \
  --species arabidopsis_thaliana --device cuda:0
python scripts/train_strand_model.py configs/publication_v1_model.json \
  --species triticum_aestivum --device cuda:0
python scripts/train_strand_model.py configs/publication_v1_model.json \
  --species zea_mays --device cuda:0
```

Activate the processing environment for the two label-execution commands;
the plan commands and pure label computations use Python's standard library.
The BAM paths are templates for the recorded processed inputs. This code release
does not assert that the Bowtie2 wrapper alone reconstructs every production BAM;
in particular, the archived GSE181488 preprocessing used HISAT2, as documented in
[production provenance](production-provenance.md).

**Preparation order matters for the archived background sample.** The code uses
`42 + species_index * 1000 + split_index` for background sampling, where species
order is Arabidopsis, wheat, maize and split order is train, validation, test.
Running `prepare --species triticum_aestivum` alone resets its species index and
samples different background windows. Use the complete configuration for the
preparation command above. Training or prediction can select one species.

After training, run genome scanning against the corresponding checkpoint:

```bash
python scripts/predict_genome.py configs/publication_v1_model.json \
  --species arabidopsis_thaliana --device cuda:0
```

The scan uses 1,024-bp inputs at 128-bp steps. Each pair of window scores is
assigned to its central 128 bp, offsets `[448, 576)` within the input. This is a
display assignment of a window prediction, not a base-resolution prediction.
Raw nonnegative plus/minus scores are retained. `minus_signed.bedGraph` negates
the minus score only for optional mirrored display. No prediction-score cutoff
or min-max normalization is applied to the saved quantitative predictions.

## Available inputs and evaluation meaning

This repository provides executable code, configuration templates, the small
coordinate masks, documented accessions/identities and the packaged figure
inputs described in the publication index. Full reference FASTA files, raw reads, processed BAMs and complete training sequence
caches are not bundled. Pretrained inference checkpoints are now available as
[release assets](https://github.com/wangyiduo217-beep/plant-pausing-bowtie2/releases/tag/strand-seiplant-v1-pretrained-20261112).
The [CPU quick start](pretrained-quickstart.md) uses bundled real DNA examples.
Original checkpoint hashes identify the archived production files; the download
manifest separately records the sanitized inference-file hashes. All model
tensors are unchanged, including the fixed B-spline buffer. Label and
manifest hashes similarly identify archived files. Fresh gzip files may have
different container metadata, even when their decompressed rows agree.

Training/evaluation metrics compare the **consensus-derived training target**
to the model score. Current signal-association figures instead pair the same
held-out 1,024-bp window and strand using `log2(1 + mean SRX CPM)` and
`log2(1 + predicted score)`. Each SRX is normalized independently before its CPM
is included in the equal-weight mean. These are separate analyses and their
correlation coefficients must not be interchanged. The older
`plot_bam_vs_prediction.py` coverage display is not the current endpoint-based
publication figure workflow.

## Upstream attribution

The architecture adapts [SeiPlant](https://github.com/compbioNJU/SeiPlant), consulted
at commit `1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853`, for the two GRO-seq regression
targets described here. Cite Lv et al., *Cross-species prediction of histone
modifications in plants via deep learning*, Genome Biology 27, 20 (2026),
[doi:10.1186/s13059-025-03929-4](https://doi.org/10.1186/s13059-025-03929-4).
Upstream SeiPlant is MIT-licensed; retain its copyright and permission notice
for incorporated/adapted portions. This model specification does not assign a
new license to the entire project.
