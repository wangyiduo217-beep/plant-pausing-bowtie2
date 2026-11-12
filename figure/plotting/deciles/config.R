# Editable appearance; scientific grouping and saved statistics are fixed.
FONT_FAMILY <- Sys.getenv("GROSEQ_FONT_FAMILY", "Arial")
FONT_PT <- 7
PANEL_PT <- 8
PNG_DPI <- 600
Y_LIMITS <- c(-0.25, 0.50)
PLUS_COLOR <- "#21844B"
MINUS_COLOR <- "#3278A5"
STRAND_OFFSET <- 0.18
POINT_SIZE_MM <- 1.25
CI_LINE_WIDTH <- 0.85
DECILE_LAYOUT <- Sys.getenv("GROSEQ_DECILE_LAYOUT", "grouped") # grouped or sixpanel
