# Open this file in RStudio and click Source, or run it with Rscript.
# This script redraws supplied Pearson r and bootstrap confidence intervals.
# It does not merge strands, regroup windows, or rerun bootstrap replicates.
# Only standard R graphics (grid and grDevices) are used.

FONT_FAMILY <- "Arial"
FONT_PT <- 7
PANEL_PT <- 8
PNG_DPI <- 600
Y_LIMITS <- c(-0.25, 0.50)
PLUS_COLOR <- "#21844B"
MINUS_COLOR <- "#3278A5"
STRAND_OFFSET <- 0.18
POINT_SIZE_MM <- 1.25
CI_LINE_WIDTH <- 0.85

if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  for (locale_name in c("English_United States.utf8", "Chinese_China.utf8")) {
    locale_result <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", locale_name), error=function(e) ""))
    if (nzchar(locale_result)) break
  }
}

find_result_dir <- function() {
  override <- Sys.getenv("GROSEQ_GROUPED_STRANDS_DIR", "")
  candidates <- if (nzchar(override)) override else character()
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f) && length(f) == 1L && file.exists(f)) {
      candidates <- c(candidates, dirname(normalizePath(f, winslash="/", mustWork=TRUE)))
    }
  }
  args <- commandArgs(trailingOnly=FALSE)
  script_args <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(script_args)) {
    candidates <- c(candidates, dirname(normalizePath(script_args[[1]], winslash="/", mustWork=TRUE)))
  }
  for (d in unique(c(candidates, getwd()))) {
    if (file.exists(file.path(d, "data", "group_statistics.csv"))) {
      return(normalizePath(d, winslash="/", mustWork=TRUE))
    }
  }
  stop("Cannot locate data/group_statistics.csv beside this script. Source the full script path, or set GROSEQ_GROUPED_STRANDS_DIR.")
}

result_dir <- find_result_dir()
if (!capabilities("cairo")) stop("Cairo graphics support is required to embed Arial text.")
statistics_path <- file.path(result_dir, "data", "group_statistics.csv")
dat <- read.csv(statistics_path, stringsAsFactors=FALSE, check.names=FALSE,
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

output_dir <- file.path(result_dir, "generated_R")
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
point_mm <- function(x, y, color) {
  grid::grid.points(x=mm(x), y=mm(y), pch=21, size=mm(POINT_SIZE_MM),
                    gp=grid::gpar(col="white", fill=color, lwd=0.35))
}

draw_panel <- function(key, column) {
  col_x <- 13.5 + (column-1)*57.5
  plot_left <- col_x + 15
  plot_width <- 47.5
  plot_bottom <- 194
  plot_height <- 44
  plot_top <- plot_bottom + plot_height
  x_position <- function(v) plot_left + (v-0.45)/(10.55-0.45)*plot_width
  x_categories <- x_position(1:10)
  y_position <- function(v) plot_bottom + (v-Y_LIMITS[1])/diff(Y_LIMITS)*plot_height

  text_mm(letters[column], plot_left-6, 264, bold=TRUE, size=PANEL_PT)
  text_mm(scientific_names[column], plot_left+plot_width/2, 264, just="center", italic=TRUE)
  line_mm(c(plot_left, plot_left), c(plot_bottom, plot_top))
  line_mm(c(plot_left, plot_left+plot_width), c(plot_bottom, plot_bottom))
  line_mm(c(plot_left, plot_left+plot_width), rep(y_position(0),2), color="#A7A7A7", width=0.55, type=2)
  for (v in seq(-0.2, 0.4, by=0.2)) {
    label <- if (abs(v) < 1e-9) "0" else format(round(v,1), nsmall=1, trim=TRUE)
    line_mm(c(plot_left-0.7, plot_left), rep(y_position(v),2))
    text_mm(label, plot_left-1.2, y_position(v), just="right")
  }
  if (column == 1L) text_mm("Pearson r", plot_left-10, plot_bottom+plot_height/2, just="center", rot=90)
  for (j in 1:10) text_mm(group_labels[j], x_categories[j], plot_bottom-2.8, just="center")

  for (strand in c("plus", "minus")) {
    d <- dat[dat$species == key & dat$strand == strand, , drop=FALSE]
    p <- d[match(1:10, d$group_id), , drop=FALSE]
    zero_n <- d$n_windows[d$group_id == 0L]
    color <- if (strand == "plus") PLUS_COLOR else MINUS_COLOR
    offset <- if (strand == "plus") -STRAND_OFFSET else STRAND_OFFSET
    count_y <- if (strand == "plus") 249 else 241
    zero_y <- if (strand == "plus") 184 else 179
    text_mm(if (strand == "plus") "n (+)" else "n (-)", plot_left-1, count_y+3.0, just="right", color=color)
    for (j in 1:10) {
      xx <- x_position(j+offset)
      lo <- y_position(p$ci95_low[j]); hi <- y_position(p$ci95_high[j])
      line_mm(rep(xx,2), c(lo,hi), color=color, width=CI_LINE_WIDTH)
      line_mm(xx+c(-0.34,0.34), c(lo,lo), color=color, width=CI_LINE_WIDTH)
      line_mm(xx+c(-0.34,0.34), c(hi,hi), color=color, width=CI_LINE_WIDTH)
      point_mm(xx, y_position(p$pearson_r[j]), color)
      text_mm(format(p$n_windows[j], big.mark=",", trim=TRUE, scientific=FALSE),
              x_categories[j], count_y, just="left", color=color, rot=90)
    }
    text_mm(paste0(if (strand == "plus") "Zero CPM (+): n = " else "Zero CPM (-): n = ",
                   format(zero_n, big.mark=",", trim=TRUE), "; r = NA"),
            plot_left+plot_width/2, zero_y, just="center", color=color)
  }
}

draw_page <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x=mm(0), y=mm(0), width=mm(page_mm[1]), height=mm(page_mm[2]),
                                    just=c("left", "bottom"), clip="off"))
  on.exit(grid::popViewport(), add=TRUE)
  text_mm("Prediction-experiment association across nonzero CPM deciles", 28.5, 282, bold=TRUE)
  text_mm("strand_seiplant_v1 | grouped strand display | statistics unchanged", 28.5, 277)
  point_mm(72, 257, PLUS_COLOR)
  text_mm("Plus strand (+)", 74.5, 257, color=PLUS_COLOR)
  point_mm(113, 257, MINUS_COLOR)
  text_mm("Minus strand (-)", 115.5, 257, color=MINUS_COLOR)
  for (column in 1:3) draw_panel(species_keys[column], column)
  text_mm("Nonzero CPM decile (Q1 lowest; Q10 highest)", 110, 172, just="center")

  text_mm("Grouped display of independently calculated strand-specific correlations.", 13.5, 155, bold=TRUE)
  caption_lines <- c(
    "a-c, Arabidopsis, wheat and maize. Within each species, each Q1-Q10 group displays separate plus- and minus-strand",
    "Pearson estimates and 95% confidence intervals. Plus-strand points are offset to the left; minus-strand points",
    "are offset to the right. The original strand-specific decile memberships, Pearson coefficients, intervals and",
    "window counts are unchanged. Positive-CPM windows were previously sorted into ten equal-count rank groups in",
    "each species and strand independently; Q1-Q10 do not imply identical absolute CPM cutoffs between strands.",
    "Both experimental mean SRX CPM and the original model score use log2(1 + value), matching the overall correlation",
    "analysis. Intervals are from the existing 1,000 chromosome-stratified paired 100-kb block-bootstrap replicates.",
    "The upper green and lower blue number rows show plus- and minus-strand window counts. Zero-CPM windows are",
    "reported separately; their correlations are undefined (NA). No strands or datasets were pooled. Arabidopsis /",
    "wheat / maize use 14 / 5 / 3 SRX, with repeat sequences retained; wheat is the original five-library model.",
    "This R script only redraws supplied statistics. It does not recalculate correlations, ranks or bootstrap."
  )
  for (i in seq_along(caption_lines)) text_mm(caption_lines[i], 13.5, 150.5-(i-1)*3.5)
}

pdf_path <- file.path(output_dir, "GROseq_CPM_deciles_grouped_strands_A4.pdf")
png_path <- file.path(output_dir, "GROseq_CPM_deciles_grouped_strands_A4.png")
grDevices::cairo_pdf(filename=pdf_path, width=page_mm[1]/25.4, height=page_mm[2]/25.4,
                     family=FONT_FAMILY, pointsize=FONT_PT, bg="white", fallback_resolution=600)
tryCatch(draw_page(), finally=grDevices::dev.off())
grDevices::png(filename=png_path, width=page_mm[1]/25.4, height=page_mm[2]/25.4,
               units="in", res=PNG_DPI, pointsize=FONT_PT, type="cairo", bg="white")
tryCatch(draw_page(), finally=grDevices::dev.off())
message("Finished: ", normalizePath(pdf_path, winslash="/", mustWork=TRUE))
message("Finished: ", normalizePath(png_path, winslash="/", mustWork=TRUE))
