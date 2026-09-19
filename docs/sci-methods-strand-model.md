# Methods: sequence-based prediction of strand-resolved nascent transcription

## Training target construction

Strand-resolved training targets were derived from the consensus GRO-seq intervals
described in [GRO-seq y-label generation](groseq-y-label-generation.md). Briefly,
independent libraries from the same species were integrated separately for the plus and
minus strands, and intervals supported by at least two independent SRX libraries were
retained. Each retained interval was assigned a reproducibility confidence between 0.1
and 1.0. Reference genomes were divided into 1,024-bp windows at 512-bp intervals. For
each strand, the continuous target was calculated as the sum, over all overlapping
strand-matched consensus intervals, of the fractional window overlap multiplied by the
interval confidence. Scores were rounded to one decimal place, clipped at 1.0, and values
below 0.1 were set to zero. The two resulting targets, denoted \(y_+\) and \(y_-\), were
retained as separate response variables.

The released target tables contained windows in which at least one strand had a non-zero
score. To expose the model to genomic background during genome-wide prediction, zero
windows were sampled from the same 1,024-bp/512-bp reference grid. Background windows
did not coincide with a positive target coordinate, did not overlap the rRNA or organellar
mask, and met the configured ambiguous-base threshold. A background-to-positive ratio of
1:1 was used independently within the training, validation, and test chromosome sets.
Sampling and record shuffling used a fixed pseudorandom seed of 42.

To prevent reference-genome input/output from limiting GPU utilization, retained windows
were ordered by genomic coordinate and encoded once into unsigned-byte sequence arrays
(`A=0`, `T=1`, `C=2`, `G=3`; one byte per base). These arrays were memory mapped during
training, and minibatch indices were shuffled with a PyTorch generator initialized with
seed 42. The sequence cache and its manifest were linked by row order and independently
checksummed.

## Sequence encoding and data partitioning

For every target window, the corresponding 1,024-bp sequence was extracted from the
forward reference strand. Bases were one-hot encoded in A, T, C, and G channel order;
ambiguous positions were represented by four zeros. Data were partitioned by chromosome
before model fitting. For *Arabidopsis thaliana*, chromosome 3 was held out for testing,
chromosome 4 for validation, and chromosomes 1, 2, and 5 for training. For *Zea mays*,
chromosome 6 was used for testing, chromosome 8 for validation, and the other assembled
nuclear chromosomes for training. These partitions follow Lv et al. (2026). For hexaploid
*Triticum aestivum*, the partition was extended by homoeologous group: chromosomes
3A/3B/3D were held out for testing, 4A/4B/4D for validation, and homoeologous groups 1,
2, 5, 6, and 7 for training. Grouping A-, B-, and D-subgenome homoeologs in the same
partition reduced leakage from highly similar homoeologous sequence. Organellar,
unplaced, and scaffold sequences were excluded.

## Neural network architecture

A separate sequence model was trained for each species. The network was adapted from
the SeiPlant architecture of Lv et al. (2026), which itself derives from Sei. The first
feature-extraction stage contained two one-dimensional convolutions with 480 channels,
kernel size 9, and padding 4, followed by a residual block of two 480-channel convolutions
with rectified linear unit (ReLU) activation. The second stage applied max pooling with
kernel and stride 4, dropout at 0.2, two convolutions that increased and maintained the
channel dimension at 640, and a two-convolution residual block. The third stage repeated
this organization with 960 channels. Five residual dilated convolutions then used 960
channels, kernel size 5, dropout at 0.1, and dilation rates 2, 4, 8, 16, and 25.

The temporal representation was projected onto 16 cubic B-spline basis functions and
flattened. Consistent with the published SeiPlant software, the spline basis was fixed and
the surrounding neural-network parameters were learned. Because the present task had only
two targets, a 256-unit ReLU layer was used instead of setting the hidden width equal to
the number of targets. A final two-unit sigmoid layer produced \(\hat y_+\) and
\(\hat y_-\). Thus, the two strands shared the sequence feature extractor but retained
separate output weights and predictions.

## Model optimization and evaluation

Models were optimized with Adam at a learning rate of \(10^{-5}\) using mean squared
error over the two continuous strand targets. Minibatches contained 256 sequences, and
automatic mixed precision was used on CUDA devices. Training continued for at most 30
epochs and stopped after five consecutive epochs without improvement in validation MSE.
The checkpoint with the lowest validation MSE was retained. No test-chromosome examples
were used for model selection. Performance on the held-out validation and test
chromosomes was summarized separately for each strand using MSE, mean absolute error,
Pearson correlation, and Spearman rank correlation.

## Genome-wide signal prediction and visualization

Each trained model was applied to the corresponding reference genome using 1,024-bp
windows advanced in 128-bp increments. To reduce sequence-boundary effects, the scalar
prediction for each input was assigned only to the central 128 bp (input positions
448–576), following the genome-track procedure of Lv et al. Consecutive predictions
therefore tiled the internal portion of each chromosome without gaps. Raw \(\hat y_+\)
and \(\hat y_-\) values were written to separate BedGraph files. Predictions were not
thresholded or min-max normalized for quantitative analysis. For mirrored genome-track
visualization only, minus-strand scores were multiplied by -1; plus-strand predictions
were plotted above zero and minus-strand predictions below zero. All statistical analyses
used the original non-negative strand scores.

## Software and reproducibility

The implementation is available in this repository under
`src/plant_pausing_bowtie2/strand_model.py`, with command-line wrappers for manifest
preparation, training, genome-wide prediction, and track plotting. Configuration files
recorded the reference genome, target release, chromosome partition, random seed, model
hyperparameters, and output directory. Generated metadata recorded input and manifest
SHA-256 hashes, the PyTorch version, best epoch, validation criterion, test metrics, and
checkpoint hash. The architecture was based on Lv T, Han Q, Li Y, et al.
*Genome Biology* 27, 20 (2026), DOI 10.1186/s13059-025-03929-4, and the public SeiPlant
implementation at commit `1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853`.
