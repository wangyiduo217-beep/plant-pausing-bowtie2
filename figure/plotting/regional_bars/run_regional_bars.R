# Run this script alone to redraw only the regional grouped bars.
# It is also called automatically from run_plots.R when MAKE_REGION_BARS is TRUE.
# No CRAN packages are needed; every value is recalculated from the supplied real windows.
run_regional_plotting <- function() {
  rb_find_dir <- function() {
    override <- Sys.getenv("GROSEQ_BUNDLE_DIR","")
    if(nzchar(override)) return(normalizePath(override,winslash="/",mustWork=TRUE))
    for(fr in rev(sys.frames())) if(!is.null(fr$ofile) && file.exists(fr$ofile)) {
      d <- dirname(normalizePath(fr$ofile,winslash="/",mustWork=TRUE))
      if(file.exists(file.path(d,"config.R"))) return(d)
    }
    f <- sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])
    if(length(f)) return(dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE)))
    if(file.exists(file.path(getwd(),"config.R"))) return(getwd())
    stop("Open the supplied R project or source the full path to plot_regional_bars.R.")
  }
  rb_dir <- rb_find_dir()
  figure_dir <- normalizePath(file.path(rb_dir,"..",".."),winslash="/",mustWork=TRUE)
  data_dir <- file.path(figure_dir,"data","regional")
  rb_config <- Sys.getenv("GROSEQ_CONFIG",file.path(rb_dir,"config.R"))
  sys.source(rb_config,envir=environment())
  stopifnot(REGION_METRIC %in% c("pearson","spearman"))
  rb_keys <- c("arabidopsis","wheat","maize")
  rb_scientific <- c("Arabidopsis thaliana","Triticum aestivum","Zea mays")
  rb_regions <- c("upstream_1_3kb","upstream_lt1kb","five_prime_utr","CDS","intron",
                  "three_prime_utr","downstream_le1kb","distal")
  rb_labels <- c("Upstream (1-3 kb)","Upstream (<1 kb)","5' UTR","CDS","Intron",
                 "3' UTR","Downstream (<=1 kb)","Distal")
  files <- file.path(data_dir,paste0(rb_keys,".regional_plot_data.tsv.gz"))
  if(!all(file.exists(files))) {
    message("Regional annotation files are not present; the other requested plots were generated.")
    return(invisible(NULL))
  }
  rb_output <- file.path(figure_dir,"results","regional_bars",paste0("regional_",REGION_METRIC,"_",Y_TRANSFORM))
  dir.create(rb_output,recursive=TRUE,showWarnings=FALSE)
  rb_transform <- function(z) {
    if(Y_TRANSFORM=="log2p1") return(log2(1+z))
    if(Y_TRANSFORM=="log2") {stopifnot(all(z>0));return(log2(z))}
    z
  }
  rb_results <- list()
  for(si in seq_along(rb_keys)) {
    key <- rb_keys[si]
    d <- read.delim(gzfile(files[si]),stringsAsFactors=FALSE)
    a <- d[,c("window_index","chrom","start0","end0","region")]
    stopifnot(nrow(d)==nrow(a),!anyDuplicated(a$window_index))
    a <- a[match(d$window_index,a$window_index),,drop=FALSE]
    stopifnot(identical(d$window_index,a$window_index),
              identical(as.character(d$chrom),as.character(a$chrom)),
              all(d$start0==a$start0),all(d$end0==a$end0),
              all(a$region %in% rb_regions))
    for(region in rb_regions) for(strand in c("plus","minus")) {
      take <- a$region==region
      x <- log2(1+d[[paste0("mean_cpm_",strand)]][take])
      y <- rb_transform(d[[paste0("prediction_",strand)]][take])
      stopifnot(all(is.finite(x)),all(is.finite(y)))
      n <- length(x)
      status <- if(n<3)"insufficient_n" else if(sd(x)==0 || sd(y)==0)"zero_variance" else "ok"
      rb_results[[length(rb_results)+1L]] <- data.frame(species=key,region=region,strand=strand,n=n,
        pearson_r=if(status=="ok")cor(x,y,method="pearson") else NA_real_,
        spearman_rho=if(status=="ok")cor(x,y,method="spearman") else NA_real_,
        status=status,y_transform=Y_TRANSFORM,stringsAsFactors=FALSE)
    }
  }
  rb_table <- do.call(rbind,rb_results)
  expected_file <- file.path(data_dir,"expected_regional_statistics_log2p1.tsv")
  if(Y_TRANSFORM=="log2p1" && file.exists(expected_file)) {
    expected <- read.delim(expected_file,stringsAsFactors=FALSE)
    keys <- function(z)paste(z$species,z$region,z$strand,sep=".")
    expected <- expected[match(keys(rb_table),keys(expected)),,drop=FALSE]
    stopifnot(all(expected$n_windows==rb_table$n),all(abs(expected$pearson_r-rb_table$pearson_r)<1e-10),
              all(abs(expected$spearman_rho-rb_table$spearman_rho)<1e-10))
  }
  write.csv(rb_table,file.path(rb_output,"regional_statistics_recalculated.csv"),row.names=FALSE,na="NA")
  rb_y_field <- if(REGION_METRIC=="pearson")"pearson_r" else "spearman_rho"
  rb_y_label <- if(REGION_METRIC=="pearson")"Pearson r" else "Spearman \u03c1"
  rb_ylim <- c(min(0,floor(min(rb_table[[rb_y_field]],na.rm=TRUE)*10)/10),1.15)
  rb_caption <- paste(
    "Grouped bars show regional correlations between experimental GRO-seq signal and sequence-model scores.",
    "Plus (green) and minus (blue) strand scores are evaluated separately on the same region-assigned windows.",
    paste0("Metric: ",REGION_METRIC,". x = log2(1 + mean SRX CPM); y transformation = ",Y_TRANSFORM,"."),
    "Each value is recalculated from the supplied original formal test-window data and per-window annotations.",
    "n labels give the number of windows in each region; they are not biological replicate counts.",
    "No invented replicate points, error bars, or confidence intervals are added.",
    "NA denotes fewer than three windows or zero variance, not a correlation of zero.",
    "Regions use the same reference-strand-aware canonical protein-coding annotation for both prediction strands.",
    "Distal is relative to this protein-coding reference. Exact definitions and annotation provenance are supplied with the data.",
    "strand_seiplant_v1; repeats retained; Arabidopsis / wheat / maize: 14 / 5 / 3 SRX.")
  writeLines(rb_caption,file.path(rb_output,"figure_legend_EN.txt"),useBytes=TRUE)
  writeLines(readLines(rb_config,warn=FALSE),file.path(rb_output,"settings_used.R"),useBytes=TRUE)
  rb_draw <- function(a4=FALSE) {
    pw <- if(a4)210 else 183;ph <- if(a4)297 else 170
    ox <- if(a4)13.5 else 0;oy <- if(a4)106 else 0
    txt <- function(x,y,z,size=FONT_PT,bold=FALSE,italic=FALSE,just="left",rot=0,col="#111111") {
      par(xpd=NA)
      grid::pushViewport(grid::viewport(clip="off"))
      grid::grid.text(z,x=grid::unit(x+ox,"mm"),y=grid::unit(y+oy,"mm"),just=just,rot=rot,
        gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=size,fontface=if(bold)"bold" else if(italic)"italic" else "plain",col=col))
      grid::popViewport()
    }
    txt(30,167,paste0("Regional ",rb_y_label," correlations"),bold=TRUE)
    cols <- c("#21844B","#3278A5");fills <- adjustcolor(cols,alpha.f=0.35)
    for(si in seq_along(rb_keys)) {
      top <- 159-(si-1)*50
      txt(1,top,letters[si],size=PANEL_PT,bold=TRUE)
      txt(30,top,rb_scientific[si],italic=TRUE)
      low <- top-26
      par(fig=c((30+ox)/pw,(181+ox)/pw,(low+oy)/ph,(low+24+oy)/ph),mai=c(0,0,0,0),
          new=TRUE,ps=FONT_PT,cex=1,cex.axis=1,las=1,xaxs="i",yaxs="i",xpd=FALSE,family=FONT_FAMILY)
      plot.new();plot.window(xlim=c(0.4,8.6),ylim=rb_ylim,xaxs="i",yaxs="i")
      axis(2,at=c(if(rb_ylim[1]<0)rb_ylim[1],0,0.5,1),tck=-0.02,mgp=c(0,0.35,0),lwd=0.5,lwd.ticks=0.5)
      segments(0.4,0,8.6,0,col="#777777",lwd=0.5)
      for(i in seq_along(rb_regions)) {
        counts <- numeric(2)
        for(j in 1:2) {
          row <- rb_table[rb_table$species==rb_keys[si] & rb_table$region==rb_regions[i] &
                            rb_table$strand==c("plus","minus")[j],,drop=FALSE]
          z <- row[[rb_y_field]];xx <- i+c(-0.18,0.18)[j];counts[j] <- row$n
          if(is.finite(z)) rect(xx-0.15,0,xx+0.15,z,col=fills[j],border=cols[j],lwd=0.7)
          else text(xx,0.04,"NA",cex=1,col=cols[j],srt=90)
        }
        stopifnot(counts[1]==counts[2])
        text(i,1.08,paste0("n=",format(counts[1],big.mark=",",trim=TRUE)),cex=1)
        xxmm <- 30+(i-0.4)/8.2*151
        txt(xxmm,low-2,rb_labels[i],just="right",rot=45)
      }
      txt(12,low+12,rb_y_label,just="center",rot=90)
    }
    for(j in 1:2) {
      xx <- 30+(j-1)*43
      grid::grid.rect(x=grid::unit(xx+1.5+ox,"mm"),y=grid::unit(9+oy,"mm"),width=grid::unit(3,"mm"),
        height=grid::unit(2,"mm"),gp=grid::gpar(fill=fills[j],col=cols[j],lwd=0.7))
      txt(xx+4.5,9,c("Plus strand","Minus strand")[j])
    }
    if(a4) {
      grid::grid.text("Correlation across genomic regions",x=grid::unit(13.5,"mm"),y=grid::unit(96,"mm"),just="left",
        gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT,fontface="bold"))
      lines <- strwrap(rb_caption,width=151)
      for(i in seq_along(lines))grid::grid.text(lines[i],x=grid::unit(13.5,"mm"),y=grid::unit(92-(i-1)*3.5,"mm"),just="left",
        gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT))
    }
  }
  for(a4 in c(FALSE,TRUE)) for(ext in c("pdf","png")) {
    target <- file.path(rb_output,paste0("regional_correlations_",if(a4)"A4" else "183x170mm",".",ext))
    ww <- if(a4)210 else 183;hh <- if(a4)297 else 170
    if(ext=="pdf")cairo_pdf(target,width=ww/25.4,height=hh/25.4,family=FONT_FAMILY,pointsize=FONT_PT,onefile=TRUE)
    else png(target,width=ww/25.4,height=hh/25.4,units="in",res=PNG_DPI,type="cairo",family=FONT_FAMILY,pointsize=FONT_PT)
    par(ps=FONT_PT,cex=1,family=FONT_FAMILY);plot.new();rb_draw(a4);dev.off()
    message("Created: ",target)
  }
  invisible(rb_table)
}
run_regional_plotting()
