# Figure 2: frozen inputs for the radar/track composite

These eight files are byte-identical to the inputs of the selected
`GROseq_radar_tracks_with_legend_A4.pdf` (2026-11-16).
Together they occupy 213,091 bytes. Their identities and the original PDF checksum
are recorded in [data_manifest.json](data_manifest.json).

[Drawing commands and panel mapping](../../plotting/figure2/README.md)

| File | Contents | Used for |
| --- | --- | --- |
| `correlation.csv` | 9 target-species/source-model pairs; pooled Pearson and sample counts | Radar a-c |
| `auroc.csv` | Strand-specific and micro summaries; native HiLearn and exact AUROC | Radar a-c use only `strand=micro`, `HiLearn_AUROC` |
| `peak_overlap.csv` | P85/P90/P95 intersection, union and confusion counts | Radar a-c use only `percentile=85`, `micro_jaccard` |
| `track_selected.csv` | Three selected regions, strand, model order, correlations and selection diagnostics | Track d-f labels and source selection |
| `track_y_limits.csv` | Fixed experimental y-axis maxima and archived track scales | Track d-f axes |
| `regions/arabidopsis_3_003700000_003800000.tracks.csv.gz` | Original observed signals, validity and model scores | Panel d |
| `regions/maize_6_119800000_119900000.tracks.csv.gz` | Original observed signals, validity and model scores | Panel e |
| `regions/wheat_3A_675800000_675900000.tracks.csv.gz` | Original observed signals, validity and model scores | Panel f |

## Quantities and coordinates

- Radar values retain their original scales. No min-max normalization, area
  scoring or averaging of the three metrics is performed. The displayed 0-1
  range covers these positive Pearson values; Pearson's theoretical range is -1 to 1.
- Global Pearson pools strand-window observations after `log2(1+x)` on the
  experimental mean CPM and predicted window score. The observation counts are
  90,132 (Arabidopsis), 99,572 (maize), and 197,312 (wheat).
- Micro AUROC pools both strands and uses `HiLearn_AUROC`, not `full_exact_AUROC`.
  P85 micro Jaccard is `TP / (TP + FP + FN)`, pooling counts; the two strands
  are thresholded independently using each signal's own percentile.
- The archived overlap table also retains P90/P95. Selection diagnostics in
  `track_selected.csv` include earlier P90 calculations. Those values do not
  replace the P85 values used for the radar axes.
- Genomic coordinates are 0-based, end-exclusive. Each input spans 1,024 bp;
  the displayed central span is `[input_start0+448, input_start0+576)` (128 bp).
- Experimental tracks display `log2(1 + mean_cpm_128_<strand>)`; prediction tracks
  show raw `prediction_<model>_<strand>` scores on 0-1 axes. The scores are not CPM.
- Regional correlations use `log2(1 + mean_cpm_1024_<strand>)` against
  `log2(1 + prediction_<model>_<strand>)` over common valid 1,024-bp inputs.
  Thus the correlation inputs differ from the experimental track's 128-bp display.
- `valid` marks usable DNA inputs. Invalid inputs remain gaps. Valid counts are
  774, 774 and 708 for the three selected regions. Own-model r is the same-species
  model's regional correlation; the two prediction tracks show cross-species models.
- Regions were selected conditional on local experimental signal and performance.
  They are illustrative examples, not an unbiased estimate of genome-wide accuracy.

This package supplies the data necessary to redraw the selected figure. Global
metric tables are frozen summaries; the eight files are not the complete raw
read data or all species-by-model per-window predictions. The original model and
held-out evaluation definitions remain documented in the figure caption and
the repository's [figure methods](../../docs/publication-figures.md).

Project-generated tables use CC BY 4.0 with attribution to GROPlant contributors.
The original line endings are preserved in Git so SHA-256 verification works
after a fresh clone or ZIP download.
