FONT_FAMILY <- Sys.getenv("GROSEQ_FONT_FAMILY", "Arial")
FONT_PT <- 7
PANEL_PT <- 8
PNG_DPI <- 600
MAKE_MAIN <- TRUE
MAKE_A4 <- TRUE
# Display windows only; all full data still enter KDE and statistical summaries.
EXPERIMENT_Y_MAX <- 4
PREDICTION_Y_MAX <- 1
# Full KDE support, deliberately independent of display windows.
EXPERIMENT_KDE_MAX <- 12
PREDICTION_KDE_MAX <- 1
VIOLIN_HALF_WIDTH <- 0.145
BANDWIDTH_ADJUST <- 1
EXPERIMENT_COLOR <- "#21844B"
PREDICTION_COLOR <- "#D55E00"
# No prediction-score threshold, min-max scaling, or within-region normalization.
