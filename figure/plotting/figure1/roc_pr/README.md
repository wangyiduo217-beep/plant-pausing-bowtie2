# Figure 1: standalone HiLearn ROC and PR reconstruction

This entry redraws the three-species ROC/PR figure from all original held-out
window labels and scores in [`figure/data/figure1/roc_inputs/`](../../../data/figure1/roc_inputs/). It calls the
unchanged HiLearn 0.2.2 curve implementation and checks its source checksum.
It also exports exact all-threshold AUROC and average precision as separate
statistics; native HiLearn PR AUC is not average precision.

From the repository root:

```sh
python -m pip install -r figure/requirements.txt
python figure/plotting/figure1/roc_pr/plot_publication_roc_pr.py --font-family "Arial"
```

Use another installed font with `--font-family`, or supply a locally licensed
TTF with `--font-file`. The default output is `figure/results/figure1/roc_pr/` and
includes figures, numerical curves and validation files. `--data-dir` and
`--output-dir` override these paths. Default paths resolve from the script,
so an absolute script path also works from another working directory.

This performs curve/statistic reconstruction from the saved labels and scores;
it does not run a model or use a GPU. The selected final workflow/evaluation
composite has its own [Figure 1 entry](../README.md), whose assembly
step uses the saved final curve vertices. Moving this standalone entry into
`figure/plotting/figure1/roc_pr/` changed its repository-root lookup and default
input/output locations, not its calculations.
