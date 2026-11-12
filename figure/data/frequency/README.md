# Experimental signal frequency references

These reference counts use every formal test window, including experimental
zeros: Arabidopsis 45,066; wheat 98,656; maize 49,786. The selected model is
repeat-retained `strand_seiplant_v1`, with 14 / 5 / 3 SRX and original five-library
wheat results.

`true_signal_frequency_bins.csv` contains 48 bins for each species/strand panel
(288 rows in total). Bins have width 0.25 on `log2(1 + mean SRX CPM)`, from 0 to 12.
They are left-closed/right-open, with the final right endpoint included. Each
percentage is `100 * window_count / n_total`; all bins in a panel sum to 100%.
The first bin includes both exact zeros and low positive values.
`frequency_summary.csv` separately records exact zero counts and percentages.

The R entry point `figure/plotting/frequency/run_signal_frequency.R` reads
the complete shared `../regional/*.regional_plot_data.tsv.gz` tables, recomputes
all frequencies and checks them against the archived bin reference. It writes
the independently calculated table as `figure/results/frequency/frequency_bins_R.csv`.
No duplicate window-level dataset is needed here.

The linear-percent main figure and logarithmic-y detail use exactly the same
counts and denominators. Empty bins are omitted from the log-axis drawing, but
remain zero in the output table. Log scaling does not modify percentages.

The distribution reflects the supplied test-window sampling design, not a
random whole-genome sample. The file manifest preserves reference identities;
no raw reads, BAMs or model checkpoints are needed for this plotting workflow.
