# Optional legacy BAM coverage comparison

This tool compares experimental BAM **coverage** with existing model prediction
bedGraphs for a selected genomic interval. It retains read 1 for paired libraries,
applies each library's RNA-strand orientation, normalizes within each SRX, and
averages SRXs equally. It requires the actual indexed BAMs and existing prediction
outputs referenced by the two configuration files.

```sh
python figure/plotting/bam_comparison/plot_bam_vs_prediction.py \
  model/configs/example_y_labels.json model/configs/example_strand_model.json \
  --species arabidopsis_thaliana --chrom 3 --start 9141194 --end 9207074 \
  --output figure/results/legacy/arabidopsis_bam_vs_prediction
```

Run from the repository root with the model environment installed. Outputs are
PNG/PDF, an aligned TSV and JSON provenance. The command reads existing model
predictions; it does not train or infer a model and does not use a GPU.

This coverage display is **not the endpoint-count signal calculation used by the
selected final publication figures**. For those results use
[combined_v2](../combined_v2/README.md) or [radar_tracks](../radar_tracks/README.md).
The optional legacy tool remains available because it has a distinct documented
function. Its move from `figure/scripts/` changes only the module search path;
the coverage calculation and plotting functions are unchanged.

For a model-only track display, use the installed
entry point `plant-pausing-model plot <model-config> ...`; see the
[model command guide](../../../model/docs/strand-seiplant-model.md).
