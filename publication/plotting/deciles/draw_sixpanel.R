# Internal drawing implementation. Run run_decile_correlations.R.
dat <- read.csv(statistics_path, fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE, check.names=FALSE,
                na.strings=c("", "NA", "NaN"))
required <- c("species", "strand", "group_id", "n_windows", "pearson_r", "ci95_low", "ci95_high")
if (!all(required %in% names(dat))) stop("The statistics CSV is missing required columns.")

species_keys <- c("arabidopsis", "wheat", "maize")
scientific_names <- c("Arabidopsis thaliana", "Triticum aestivum", "Zea mays")
expected_totals <- c(arabidopsis=45066, wheat=98656, maize=49786)
group_labels <- paste0("Q", 1:10)
for (key in species_keys) for (strand in c("plus", "minus")) {
  d <- dat[dat$species == key & dat$strand == strand, , drop=FALSE]
  stopifnot(nrow(d) == 11L, identical(sort(as.integer(d$group_id)), 0:10),
            sum(d$n_windows) == expected_totals[[key]])
  z <- d[d$group_id == 0L, , drop=FALSE]
  stopifnot(is.na(z$pearson_r), is.na(z$ci95_low), is.na(z$ci95_high))
  p <- d[d$group_id > 0L, , drop=FALSE]
  stopifnot(all(is.finite(p$pearson_r)), all(is.finite(p$ci95_low)),
            all(is.finite(p$ci95_high)), all(p$ci95_low <= p$ci95_high))
  if (min(p$ci95_low) < Y_LIMITS[1] || max(p$ci95_high) > Y_LIMITS[2]) {
    stop("Y_LIMITS would clip a supplied confidence interval. Increase Y_LIMITS.")
  }
}

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
page_mm <- c(210, 297)
mm <- function(x) grid::unit(x, "mm")
text_mm <- function(label, x, y, just="left", bold=FALSE, italic=FALSE,
                    size=FONT_PT, color="#222222", rot=0) {
  face <- if (bold && italic) "bold.italic" else if (bold) "bold" else if (italic) "italic" else "plain"
  grid::grid.text(label, x=mm(x), y=mm(y), just=just, rot=rot,
                  gp=grid::gpar(fontfamily=FONT_FAMILY, fontsize=size, fontface=face, col=color))
}
line_mm <- function(x, y, color="#222222", width=0.55, type=1) {
  grid::grid.lines(x=mm(x), y=mm(y), gp=grid::gpar(col=color, lwd=width, lty=type, lineend="butt"))
}

draw_panel <- function(key, strand, column, row) {
  d <- dat[dat$species == key & dat$strand == strand, , drop=FALSE]
  zero_n <- d$n_windows[d$group_id == 0L]
  p <- d[match(1:10, d$group_id), , drop=FALSE]
  col_x <- 13.5 + (column-1)*61
  plot_left <- col_x + 12
  plot_width <- 47
  plot_bottom <- if (row == 1L) 200 else 124
  plot_height <- 46
  plot_top <- plot_bottom + plot_height
  x_positions <- plot_left + ((1:10)-0.5)/10*plot_width
  y_position <- function(v) plot_bottom + (v-Y_LIMITS[1])/diff(Y_LIMITS)*plot_height
  color <- if (strand == "plus") PLUS_COLOR else MINUS_COLOR
  letter <- letters[(row-1L)*3L + column]

  text_mm(letter, col_x+6, plot_top+17.0, bold=TRUE, size=PANEL_PT)
  text_mm(if (strand == "plus") "Plus strand (+)" else "Minus strand (-)",
          plot_left+plot_width/2, plot_top+17, just="center", color=color)
  text_mm("n", plot_left-1.0, plot_top+6.0, just="right", color="#555555")
  for (j in 1:10) {
    text_mm(format(p$n_windows[j], big.mark=",", trim=TRUE, scientific=FALSE),
            x_positions[j], plot_top+3.0, just="left", color="#555555", rot=90)
  }

  line_mm(c(plot_left, plot_left), c(plot_bottom, plot_top))
  line_mm(c(plot_left, plot_left+plot_width), c(plot_bottom, plot_bottom))
  line_mm(c(plot_left, plot_left+plot_width), rep(y_position(0),2), color="#A7A7A7", width=0.55, type=2)
  for (v in seq(-0.2, 0.4, by=0.2)) {
    label <- if (abs(v) < 1e-9) "0" else sprintf("%.1f", v)
    line_mm(c(plot_left-0.7, plot_left), rep(y_position(v),2))
    text_mm(label, plot_left-1.2, y_position(v), just="right")
  }
  if (column == 1L) text_mm("Pearson r", col_x+3.2, plot_bottom+plot_height/2, just="center", rot=90)

  for (j in 1:10) {
    lo <- y_position(p$ci95_low[j]); hi <- y_position(p$ci95_high[j])
    line_mm(rep(x_positions[j],2), c(lo,hi), color=color, width=CI_LINE_WIDTH)
    line_mm(x_positions[j]+c(-0.8,0.8), c(lo,lo), color=color, width=CI_LINE_WIDTH)
    line_mm(x_positions[j]+c(-0.8,0.8), c(hi,hi), color=color, width=CI_LINE_WIDTH)
    grid::grid.points(x=mm(x_positions[j]), y=mm(y_position(p$pearson_r[j])),
                      pch=21, size=mm(POINT_SIZE_MM),
                      gp=grid::gpar(col="white", fill=color, lwd=0.35))
    text_mm(group_labels[j], x_positions[j], plot_bottom-2.5, just="center")
  }
  text_mm(paste0("Zero-signal n = ", format(zero_n, big.mark=",", trim=TRUE), "; r = NA"),
          plot_left+plot_width/2, plot_bottom-8.0, just="center", color="#666666")
}

draw_page <- function(crop=FALSE) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x=mm(if(crop) -13.5 else 0), y=mm(if(crop) -106 else 0), width=mm(page_mm[1]), height=mm(page_mm[2]),
                                    just=c("left", "bottom"), clip="off"))
  on.exit(grid::popViewport(), add=TRUE)
  text_mm("Prediction-experiment association across CPM deciles", 25.5, 284, bold=TRUE)
  text_mm("strand_seiplant_v1 | 1,024-bp test windows | repeat sequences retained", 25.5, 279)
  for (column in 1:3) {
    text_mm(scientific_names[column], 13.5+(column-1)*61+35.5, 269, just="center", italic=TRUE)
    draw_panel(species_keys[column], "plus", column, 1L)
    draw_panel(species_keys[column], "minus", column, 2L)
  }
  text_mm("Positive-CPM rank group (Q1 lowest; Q10 highest)", 110, 110, just="center")

  text_mm("CPM-decile association in held-out test windows.", 13.5, 99, bold=TRUE)
  caption_lines <- c(
    "a-c, Plus strand; d-f, minus strand. Columns show Arabidopsis, wheat and maize. Within each species and strand,",
    "positive-CPM windows were sorted into ten equal-count rank groups (Q1-Q10). Equal CPM values can span adjacent",
    "groups; chromosome, coordinates and original window index resolve ties deterministically. Zero-CPM windows are",
    "reported separately. Points show Pearson r across all windows in each rank group. Both experimental mean SRX",
    "CPM and the original prediction score use log2(1 + value), matching the overall correlation analysis. Whiskers",
    "are percentile 95% intervals from 1,000 paired, chromosome-stratified 100-kb block-bootstrap replicates; group",
    "memberships are fixed. The same block draws are used across groups and strands. Rotated numbers show window n.",
    "Only occupied test-window blocks are sampled. Zero-signal correlations are NA because their true values are",
    "constant. Arabidopsis / wheat / maize: 45,066 / 98,656 / 49,786 windows and 14 / 5 / 3 SRX; wheat is the original",
    "five-library model. Model scores are consensus-derived 0-1 scores, not CPM estimates. Grouping restricts the",
    "signal range, so correlations are descriptive conditional associations, not causal effects of coverage. This R",
    "script redraws supplied summary statistics; it does not recompute correlations or bootstrap replicates."
  )
  for (i in seq_along(caption_lines)) text_mm(caption_lines[i], 13.5, 94.5-(i-1)*3.5)
}


for (crop in c(FALSE, TRUE)) {
  dimensions <- if(crop) c(183,183) else c(210,297)
  stem <- paste0("CPM_decile_correlations_sixpanel",if(crop) "" else "_A4")
  pdf_path <- file.path(output_dir,paste0(stem,".pdf"))
  png_path <- file.path(output_dir,paste0(stem,".png"))
  grDevices::cairo_pdf(pdf_path,width=dimensions[1]/25.4,height=dimensions[2]/25.4,
    family=FONT_FAMILY,pointsize=FONT_PT,bg="white",fallback_resolution=600)
  tryCatch(draw_page(crop=crop),finally=grDevices::dev.off())
  grDevices::png(png_path,width=dimensions[1]/25.4,height=dimensions[2]/25.4,
    units="in",res=PNG_DPI,pointsize=FONT_PT,type="cairo",bg="white")
  tryCatch(draw_page(crop=crop),finally=grDevices::dev.off())
  message("Finished: ",normalizePath(pdf_path,winslash="/",mustWork=TRUE))
}
