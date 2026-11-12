# Current extended track examples

The three `*.selected_track.csv` files contain the frozen display bins for Arabidopsis (40 kb), wheat (80 kb) and maize (80 kb). `selected_regions.csv` gives the exact 0-based half-open interval, sample sizes, selection metadata and local track statistics for each example. Full formal-test statistics shown alongside the tracks are recalculated from `../regional/` instead.

The table schema is listed per file in `data_manifest.json`. `start0`/`end0` refer to the central 128 bp display bin, `input_start0` to the 1,024 bp DNA input, `gro_log2p1_plus/minus` to the true 128 bp-bin GRO-seq signal, and `prediction_plus/minus` to raw 0-1 model scores. The true signal is `log2(1 + mean independently normalized SRX CPM)` using the same old 14/5/3 library sets. `valid` marks whether the model input is valid. Missing prediction values must remain NA, not zero.

The selected intervals are chr3:4,610,000-4,650,000, chr3A:292,140,000-292,220,000 and chr6:74,920,000-75,000,000, using 0-based half-open coordinates. The wheat interval has 24 invalid prediction bins out of 624 bins. The original left boundaries are preserved, and wheat/maize extend farther right than their earlier 40 kb examples.

`data_manifest.json` gives filenames relative to `figure/` and SHA-256 checksums. These are actual display inputs, not synthetic example values. Predictions are window scores placed on central bins; the display is not a base-resolution model output.
