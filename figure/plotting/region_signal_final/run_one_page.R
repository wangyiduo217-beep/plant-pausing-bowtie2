# Standard R only. Source this file from any working directory.
# Data paths are relative to this script: ../../data/regional/*.regional_plot_data.tsv.gz.
find_bundle_dir <- function() {
  override <- Sys.getenv("GROSEQ_ONE_PAGE_DIR", "")
  if (nzchar(override)) return(normalizePath(override, winslash="/", mustWork=TRUE))
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile) && file.exists(fr$ofile)) {
    candidate <- dirname(normalizePath(fr$ofile, winslash="/", mustWork=TRUE))
    if (file.exists(file.path(candidate, "config.R")) &&
        file.exists(file.path(candidate, "run_one_page.R"))) return(candidate)
  }
  args <- commandArgs(FALSE)
  paths <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(paths)) return(dirname(normalizePath(paths[[1]], winslash="/", mustWork=TRUE)))
  if (file.exists("run_one_page.R") && file.exists("config.R")) return(getwd())
  stop("Source the full path to run_one_page.R or open the supplied project.")
}
bundle <- find_bundle_dir()
sys.source(file.path(bundle, "config.R"), envir=environment())
stopifnot(capabilities("cairo"), CORRELATION_Y_MAX == 1, MEAN_CPM_Y_MAX > 0)
data_dir <- normalizePath(file.path(bundle, "..", "..", "data", "regional"), winslash="/", mustWork=TRUE)
output_dir <- file.path(bundle, "..", "..", "results", "region_signal_final_R")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

species_keys <- c("arabidopsis", "wheat", "maize")
species_names <- c("Arabidopsis thaliana", "Triticum aestivum", "Zea mays")
regions <- c("upstream_1_3kb", "upstream_lt1kb", "five_prime_utr", "CDS", "intron",
             "three_prime_utr", "downstream_le1kb", "distal")
region_labels <- c("Up\n1-3kb", "Up\n<1kb", "5'UTR\n", "CDS\n", "Intron\n",
                   "3'UTR\n", "Down\n<=1kb", "Distal\n")
expected_n <- c(45066L, 98656L, 49786L)
chromosomes <- c("3", "3A / 3B / 3D", "6")
srx_counts <- c(14L, 5L, 3L)

# Each row is a complete region, species, and strand. Pearson r and Spearman
# rho use log2(1 + per-window mean SRX CPM) versus log2(1 + predicted score).
# The plotted right-axis quantity is mean(raw per-window mean CPM), including
# all exact zeros. No prediction value enters this right-axis average.
records <- list()
region_counts <- list()
for (si in seq_along(species_keys)) {
  key <- species_keys[si]
  input <- file.path(data_dir, paste0(key, ".regional_plot_data.tsv.gz"))
  d <- read.delim(gzfile(input), stringsAsFactors=FALSE)
  stopifnot(nrow(d) == expected_n[si], all(d$region %in% regions),
            all(d$window_index == seq_len(nrow(d))))
  region_counts[[key]] <- tabulate(match(d$region, regions), nbins=length(regions))
  for (strand in c("plus", "minus")) {
    raw_cpm <- d[[paste0("mean_cpm_", strand)]]
    raw_prediction <- d[[paste0("prediction_", strand)]]
    stopifnot(all(is.finite(raw_cpm)), all(raw_cpm >= 0),
              all(is.finite(raw_prediction)), all(raw_prediction >= 0), all(raw_prediction <= 1))
    experimental_log <- log2(1 + raw_cpm)
    prediction_log <- log2(1 + raw_prediction)
    for (region in regions) {
      take <- d$region == region
      n <- sum(take)
      x <- experimental_log[take]
      y <- prediction_log[take]
      stopifnot(n > 1L, sd(x) > 0, sd(y) > 0)
      records[[length(records) + 1L]] <- data.frame(
        species=key, species_name=species_names[si], strand=strand, region=region,
        n_windows=n, raw_cpm_zero_count=sum(raw_cpm[take] == 0),
        raw_cpm_mean=mean(raw_cpm[take]),
        experimental_log2p1_mean=mean(x), prediction_raw_mean=mean(raw_prediction[take]),
        pearson_r=cor(x, y, method="pearson"),
        spearman_rho=cor(x, y, method="spearman"), stringsAsFactors=FALSE)
    }
  }
}
statistics <- do.call(rbind, records)
rownames(statistics) <- NULL
stopifnot(nrow(statistics) == 48L, sum(expected_n) == 193508L,
          sum(statistics$n_windows) == 2L * sum(expected_n),
          all(is.finite(statistics$pearson_r)), all(statistics$pearson_r >= 0),
          all(statistics$pearson_r <= CORRELATION_Y_MAX),
          all(statistics$raw_cpm_mean <= MEAN_CPM_Y_MAX))
options(digits=17)
write.table(statistics, file.path(output_dir, "recalculated_region_statistics.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Optional numerical regression check against the unchanged supplied summary.
reference_file <- file.path(bundle, "..", "..", "data", "region_signal_final", "regional_correlation_and_mean_signal.tsv")
if (file.exists(reference_file)) {
  reference <- read.delim(reference_file, stringsAsFactors=FALSE)
  make_key <- function(d) paste(d$species, d$strand, d$region, sep="/")
  stopifnot(nrow(reference) == 48L, !anyDuplicated(make_key(reference)))
  reference <- reference[match(make_key(statistics), make_key(reference)), ]
  stopifnot(!anyNA(reference$region))
  matching_columns <- c(n_windows="n_windows", raw_cpm_zero_count="zero_count",
                        raw_cpm_mean="mean_cpm", experimental_log2p1_mean="mean_log2p1_cpm",
                        pearson_r="pearson_r", spearman_rho="spearman_rho")
  for (column in names(matching_columns)) {
    stopifnot(max(abs(statistics[[column]] - reference[[matching_columns[[column]]]])) <= 1e-12)
  }
}

ctx <- list(width=183, height=170, ox=0, oy=0)
PLOT_X <- c(12, 110)
PLOT_WIDTH <- 64
PLOT_HEIGHT <- 28.5
X_LIMITS <- c(-0.6, 7.6)
text_at <- function(x, y, label, size=FONT_PT, bold=FALSE, italic=FALSE,
                    just="left", col="#111111", baseline=TRUE) {
  if (baseline && length(just) == 1L && just %in% c("left", "right", "center")) {
    y <- y + size * 25.4 / 72 * 0.27
  }
  grid::pushViewport(grid::viewport(clip="off"))
  grid::grid.text(label, x=grid::unit(x + ctx$ox, "mm"),
                  y=grid::unit(y + ctx$oy, "mm"), just=just,
                  gp=grid::gpar(fontfamily=FONT_FAMILY, fontsize=size,
                    fontface=if (bold) "bold" else if (italic) "italic" else "plain", col=col))
  grid::popViewport()
}
rotated_text <- function(x, y, label, col) {
  grid::pushViewport(grid::viewport(clip="off"))
  grid::grid.text(label, x=grid::unit(x + ctx$ox, "mm"),
                  y=grid::unit(y + ctx$oy, "mm"), rot=90,
                  gp=grid::gpar(fontfamily=FONT_FAMILY, fontsize=FONT_PT, col=col))
  grid::popViewport()
}
group_x <- function(j, column) PLOT_X[column] + (j - 1 - X_LIMITS[1]) / diff(X_LIMITS) * PLOT_WIDTH
format_tick <- function(v) vapply(v, function(value) format(value, trim=TRUE, scientific=FALSE), character(1))
correlation_label <- function(x, y, label, col) {
  grid::pushViewport(grid::viewport(clip="off"))
  gp <- grid::gpar(fontfamily=FONT_FAMILY, fontsize=FONT_PT, col=col)
  grob <- grid::textGrob(label, gp=gp)
  yy <- y + FONT_PT * 25.4 / 72 * 0.27
  grid::grid.rect(x=grid::unit(x + ctx$ox, "mm"), y=grid::unit(yy + ctx$oy, "mm"),
                  width=grid::grobWidth(grob) + grid::unit(0.25, "mm"),
                  height=grid::grobHeight(grob) + grid::unit(0.15, "mm"),
                  gp=grid::gpar(fill="white", col=NA))
  grid::grid.text(label, x=grid::unit(x + ctx$ox, "mm"), y=grid::unit(yy + ctx$oy, "mm"), gp=gp)
  grid::popViewport()
}

# Inspect the complete line segment under each label, including adjacent slopes.
# Only label height changes. Correlations and raw CPM coordinates stay unchanged.
label_positions <- function(r, cpm) {
  label_y <- r + 0.025
  line_height <- cpm / MEAN_CPM_Y_MAX
  for (j in seq_along(r)) {
    xj <- j - 1
    local_x <- sort(unique(c(max(0, xj - 0.35), xj, min(length(r) - 1, xj + 0.35))))
    local_line <- approx(seq_along(r) - 1, line_height, xout=local_x)$y
    local_range <- range(local_line)
    if (local_range[2] >= label_y[j] - 0.035 && local_range[1] <= label_y[j] + 0.11) {
      label_y[j] <- max(label_y[j], local_range[2] + 0.055)
    }
  }
  stopifnot(all(label_y < 0.93))
  label_y
}

draw_panel <- function(si, strand, column, y, h=PLOT_HEIGHT) {
  d <- statistics[statistics$species == species_keys[si] & statistics$strand == strand, ]
  d <- d[match(regions, d$region), ]
  bar_color <- if (strand == "plus") PLUS_COLOR else MINUS_COLOR
  x <- PLOT_X[column]
  par(fig=c((x + ctx$ox) / ctx$width, (x + PLOT_WIDTH + ctx$ox) / ctx$width,
            (y + ctx$oy) / ctx$height, (y + h + ctx$oy) / ctx$height),
      mai=c(0, 0, 0, 0), new=TRUE, ps=FONT_PT, cex=1, cex.axis=1,
      cex.lab=1, family=FONT_FAMILY, xaxs="i", yaxs="i", xpd=FALSE)
  plot.new()
  plot.window(X_LIMITS, c(0, CORRELATION_Y_MAX), xaxs="i", yaxs="i")
  left_usr <- par("usr")
  axis(2, at=CORRELATION_TICKS, labels=c("0", "0.25", "0.50", "0.75", "1"), las=1,
       tck=-0.018, lwd=0.5, mgp=c(0, 0.3, 0), cex.axis=1,
       col=bar_color, col.axis=bar_color)
  segments(X_LIMITS[1], 0, X_LIMITS[2], 0, col="#888888", lwd=0.5)
  xpos <- seq_along(regions) - 1
  rect(xpos - BAR_WIDTH / 2, 0, xpos + BAR_WIDTH / 2, d$pearson_r,
       col=adjustcolor(bar_color, BAR_ALPHA), border=bar_color, lwd=BAR_BORDER_WIDTH)
  for (j in seq_along(regions)) {
    text_at(group_x(j, column), y - 1.4,
            paste0(region_labels[j], "\n", format(d$n_windows[j], trim=TRUE, scientific=FALSE)),
            just="top")
  }
  # The native right axis receives raw arithmetic means directly.
  par(usr=c(X_LIMITS, 0, MEAN_CPM_Y_MAX), xpd=FALSE)
  stopifnot(isTRUE(all.equal(unname(par("usr")), c(X_LIMITS, 0, MEAN_CPM_Y_MAX))))
  axis(4, at=MEAN_CPM_TICKS, labels=format_tick(MEAN_CPM_TICKS), las=1,
       tck=-0.018, lwd=0.5, mgp=c(0, 0.3, 0), cex.axis=1,
       col=MEAN_COLOR, col.axis=MEAN_COLOR)
  lines(xpos, d$raw_cpm_mean, col=MEAN_COLOR, lwd=LINE_WIDTH)
  points(xpos, d$raw_cpm_mean, col=MEAN_COLOR, pch=16, cex=POINT_CEX)
  par(usr=left_usr, xpd=FALSE)
  rotated_text(x - 10.4, y + h / 2, "Pearson r", bar_color)
  rotated_text(x + PLOT_WIDTH + 8, y + h / 2, "Mean CPM", MEAN_COLOR)
  label_y <- label_positions(d$pearson_r, d$raw_cpm_mean)
  for (j in seq_along(regions)) {
    correlation_label(group_x(j, column), y + h * label_y[j], sprintf("%.2f", d$pearson_r[j]), bar_color)
  }
}

draw_legend <- function() {
  for (j in seq_along(c(PLUS_COLOR, MINUS_COLOR))) {
    color <- c(PLUS_COLOR, MINUS_COLOR)[j]
    grid::grid.rect(x=grid::unit(114 + (j - 1) * 2.2 + ctx$ox, "mm"),
                    y=grid::unit(168.1 + ctx$oy, "mm"),
                    width=grid::unit(1.8, "mm"), height=grid::unit(1.8, "mm"),
                    gp=grid::gpar(fill=adjustcolor(color, BAR_ALPHA), col=color, lwd=BAR_BORDER_WIDTH))
  }
  text_at(118.5, 168.1, "Pearson r", baseline=FALSE)
  grid::grid.lines(x=grid::unit(c(139.5, 143) + ctx$ox, "mm"),
                   y=grid::unit(c(168.1, 168.1) + ctx$oy, "mm"),
                   gp=grid::gpar(col=MEAN_COLOR, lwd=LINE_WIDTH))
  grid::grid.points(x=grid::unit(141.25 + ctx$ox, "mm"),
                    y=grid::unit(168.1 + ctx$oy, "mm"), pch=16,
                    size=grid::unit(1.3, "mm"), gp=grid::gpar(col=MEAN_COLOR))
  text_at(145, 168.1, "Mean CPM", baseline=FALSE)
}

draw_caption <- function() {
  title <- "Regional prediction agreement and mean experimental GRO-seq signal: plus and minus strands."
  lines <- c(
    "Bars show Pearson r between log2(1 + per-window mean SRX CPM) and log2(1 + predicted score) within each region.",
    "Orange points and lines show the arithmetic mean of raw per-window mean CPM, including every zero-valued window.",
    "The right axis shows raw CPM, not mean log2 signal or predicted scores. Both axes use fixed ranges: r = 0-1; CPM = 0-6.",
    "All 193,508 original held-out windows contribute separately to each strand, without score thresholds or outlier filtering.",
    "Up = upstream; Down = downstream. Labels give distance to the gene boundary; categories use window-midpoint annotations.",
    "Counts below categories are the complete regional window counts. Values above bars are Pearson r rounded to two decimals.",
    "Categories are not a continuous genomic-distance or Mb axis. The two axes have independent units; visual correspondence",
    "is descriptive and does not establish causation. No significance tests are performed. Tables retain full numerical precision.")
  text_at(0, -10, title, bold=TRUE)
  for (i in seq_along(lines)) text_at(0, -14.1 - (i - 1) * 3.5, lines[i])
}

draw_current <- function(a4=FALSE) {
  text_at(0, 168.1, "Regional correlation and mean GRO-seq signal", bold=TRUE, baseline=FALSE)
  draw_legend()
  text_at(44, 157, "Plus strand (+)", bold=TRUE, just="center", col=PLUS_COLOR)
  text_at(142, 157, "Minus strand (-)", bold=TRUE, just="center", col=MINUS_COLOR)
  for (si in seq_along(species_keys)) {
    header <- c(149, 99, 49)[si]
    text_at(0, header, species_names[si], italic=TRUE)
    text_at(183, header, paste0("chr", chromosomes[si], "; n = ",
                              format(expected_n[si], big.mark=",", trim=TRUE), "; ", srx_counts[si], " SRX"),
            just="right")
    for (column in 1:2) {
      strand <- c("plus", "minus")[column]
      text_at(c(0, 98)[column], header - 6, letters[(si - 1) * 2 + column], size=PANEL_PT, bold=TRUE)
      draw_panel(si, strand, column, y=header - 34)
    }
  }
  text_at(0, 2.8, "Up = upstream; Down = downstream; counts = windows. strand_seiplant_v1; repeats retained; zeros included.")
  if (a4) draw_caption()
}

open_page <- function(file, a4=FALSE) {
  width <- if (a4) 210 else 183
  height <- if (a4) 297 else 170
  if (grepl("\\.pdf$", file)) {
    cairo_pdf(file, width=width / 25.4, height=height / 25.4,
              pointsize=FONT_PT, family=FONT_FAMILY, onefile=TRUE)
  } else {
    png(file, width=width / 25.4, height=height / 25.4, units="in", res=PNG_DPI,
        type="cairo", pointsize=FONT_PT, family=FONT_FAMILY)
  }
  ctx <<- list(width=width, height=height, ox=if (a4) 13.5 else 0, oy=if (a4) 106 else 0)
  par(ps=FONT_PT, cex=1, cex.axis=1, cex.lab=1, family=FONT_FAMILY)
  plot.new()
}

for (a4 in c(FALSE, TRUE)) {
  if ((a4 && !MAKE_A4) || (!a4 && !MAKE_MAIN)) next
  suffix <- if (a4) "A4" else "183x170mm"
  file <- file.path(output_dir, paste0("regional_correlation_mean_signal_one_page_", suffix, ".pdf"))
  open_page(file, a4)
  draw_current(a4)
  dev.off()
  cat("Created:", file, "\n")
}
if (MAKE_MAIN) {
  file <- file.path(output_dir, "regional_correlation_mean_signal_one_page_183x170mm.png")
  open_page(file, FALSE)
  draw_current(FALSE)
  dev.off()
  cat("Created:", file, "\n")
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
cat("Validated: all 193508 windows per strand; 48 full-region statistics; unchanged log2p1 correlations and raw CPM means including zeros; native independent axes; one page per PDF.\n")