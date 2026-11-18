# Exact inputs for the final combined figure

These files belong to `GROseq_combined_v2_A4.pdf` (2026-11-12), not to the later
cross-species radar figure. Inputs are copied byte-for-byte from the selected
original result's data sources.

- `test_windows/`: three `*.all_test_windows.npz` files, containing every held-out
  test window's original predictions, mean experimental CPM and underlying SRX
  endpoint counts/normalization denominators. Counts are 45,066 Arabidopsis,
  98,656 wheat and 49,786 maize. Both strands are retained. The accompanying
  legacy correlation/OLS TSV is a provenance check; the plotted transformed
  statistics are calculated from the arrays by `build_source_panels.py`.
- `tracks/`: three `*.selected_track.npz` files and `selected_regions.tsv` for
  the exact representative regions. Each model prediction is for a 1,024-bp
  input and displayed in its central 128-bp bin. Invalid DNA remains missing,
  with the original gray unavailable-input regions; it is not replaced by zero.
- `curves/`: three `*.HiLearn_curves.npz` files containing the exact ROC and PR
  vertices used in the selected figure. Both output strands were pooled for
  these micro curves. PR area is native HiLearn trapezoidal area, including
  its artificial endpoints, not average precision.

`source_data_manifest.json` records the eight original source-panel input hashes
relative to this directory. Curve and vector-asset hashes are recorded in
[`source_manifest.json`](../../plotting/combined_v2/provenance/source_manifest.json).
The primary redraw scripts and methods are
[`figure/plotting/combined_v2`](../../plotting/combined_v2/).
