# Change these settings, save, then run source("run_plots.R", encoding="UTF-8").
# Only base R and the standard R graphics capabilities are required.
POINT_STYLE <- "density"      # "density", "scatter", or "hexbin"
FIT_METHOD <- "ols"           # "ols", "local_linear", or "none"
Y_TRANSFORM <- "log2p1"       # "log2p1", "log2", or "raw"
MIN_PREDICTION_SCORE <- 0 # Keep raw score >= this value, separately for each strand.
COMPARISON_THRESHOLDS <- c(0, 0.001, 0.005, 0.01)
X_LIMITS <- c(0, 12.2)        # No negative axis margins in this comparison.
Y_LIMITS <- c(0, 1.02)        # log2(1 + raw score); no coordinate shift.
POINT_CEX <- 0.17
POINT_ALPHA <- 0.8
SCATTER_COLOR <- "#9a9a9a"
DENSITY_PALETTE <- c("#c4c4c4", "#929292", "#4b4b4b", "#111111")
HEX_GRIDSIZE <- c(36L, 21L)
FIT_COLOR <- "#555555"
FIT_LINE_TYPE <- 2            # Dashed straight line, as in the first version.
FIT_LINE_WIDTH <- 0.8
STAT_COLOR <- "#D55E00"
FONT_FAMILY <- Sys.getenv("GROSEQ_FONT_FAMILY", "Arial")
FONT_PT <- 7
PANEL_PT <- 8
PNG_DPI <- 600
LOCAL_BANDWIDTH <- 0.6
MIN_EFFECTIVE_WEIGHT_COUNT <- 30
MAKE_INDIVIDUAL_SPECIES <- FALSE
MAKE_TRACK_COMPOSITE <- FALSE # Optional; raw Track values are never filtered.
MAKE_A4_PROOFS <- TRUE
MAKE_REGION_BARS <- FALSE     # This small bundle focuses on threshold scatter comparisons.
REGION_METRIC <- "pearson"     # "pearson" or "spearman"

DRAW_CORRELATIONS <- TRUE
OUTPUT_SUBDIR <- "correlations"
