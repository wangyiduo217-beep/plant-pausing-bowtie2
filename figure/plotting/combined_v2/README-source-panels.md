# Data-driven source figure for the final GRO-seq combined figure

This self-contained bundle redraws the source tracks and all-held-out-window scatter panels used by `GROseq_combined_v2_A4.pdf`. It is a path/font portability adaptation of the archived original, not an image reconstruction. Original input files are copied byte-for-byte; see `../../data/combined_v2/source_data_manifest.json`.

```
python figure/plotting/combined_v2/build_source_panels.py
```

Install numpy, scipy, matplotlib, reportlab, pypdf and pdfplumber. Exact typography requires licensed Arial regular/bold/italic/bold italic installed on the system, or `--font-dir /path/to/Arial`. No font binaries are distributed. All coordinates, scores, density calculations and statistical methods are preserved. The six scatter panels draw every original test window and retain their original rasterized point layers; axes, labels and OLS lines are vector artwork.

Outputs include `figure/results/combined_v2_source/GROseq_v1_tracks_OLS_correlations_A4.pdf`, the 183 x 170 mm source artwork, the standalone six-scatter figure, exact calculated plotting arrays, and structural/numerical validation. The final combined figure recolors/repositions some annotations separately. Numerical reproduction does not imply byte-identical PDF serialization across library versions.

## Exact plot inputs

* test_windows/*.all_test_windows.npz: original held-out windows, predictions, mean CPM, and underlying per-SRX endpoint counts/normalization denominators. No subsampling.
* tracks/*.selected_track.npz: original 128-bp display bins; full 1024-bp-input predictions shown in the central 128 bp. Missing wheat DNA inputs remain NaN and gray-shaded.
* tracks/selected_regions.tsv: original displayed genomic regions.
* test_windows/correlation_and_OLS.tsv: legacy raw-score statistics used for source consistency, not reused as the transformed plotted correlation.

Per-species formal test-window counts are 45,066 Arabidopsis, 98,656 wheat, and 49,786 maize; both strands are shown. The new figure concerns the legacy strand_seiplant_v1 model with 14/5/3 SRX in that species order. It does not switch to a later wheat model.

The archived original is `provenance/original_build_final.py`. Model training, model inference, BAM processing and background attribution calculations are not run by this plotting script.

To assemble the final composite from these newly rendered source panels:

```sh
python figure/plotting/combined_v2/build_combined_v2.py --source-figure2 figure/results/combined_v2_source/GROseq_v1_tracks_OLS_correlations_A4.pdf
```
