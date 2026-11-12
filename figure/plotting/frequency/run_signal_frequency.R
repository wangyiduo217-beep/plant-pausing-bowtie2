# Independently recompute experimental histograms from complete shared window tables.
find_module_dir <- function(filename) {
  candidates <- character()
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile) && basename(fr$ofile) == filename)
    candidates <- c(candidates, dirname(normalizePath(fr$ofile,winslash="/",mustWork=TRUE)))
  args <- commandArgs(FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if(length(f) && basename(f[[1]]) == filename)
    candidates <- c(candidates, dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE)))
  for (d in unique(candidates)) if(file.exists(file.path(d, filename))) return(d)
  stop("Source the full path to ", filename)
}

module_dir <- find_module_dir("run_signal_frequency.R")
figure_dir <- normalizePath(file.path(module_dir,"../.."),winslash="/",mustWork=TRUE)
source(file.path(module_dir,"config.R"),local=environment(),encoding="UTF-8")
input_dir <- file.path(figure_dir,"data/regional")
output_dir <- file.path(figure_dir,"results/frequency")
if(!capabilities("cairo"))stop("Cairo graphics support is required for embedded text.")
if(BIN_WIDTH != 0.25 || X_MAX != 12)stop("This frozen analysis uses bin width 0.25 and x range 0-12.")
species_keys <- c("arabidopsis","wheat","maize")
scientific_names <- c("Arabidopsis thaliana","Triticum aestivum","Zea mays")
expected_n <- c(45066L,98656L,49786L)
edges <- seq(0,X_MAX,by=BIN_WIDTH)
panels <- list()
bin_rows <- list()
counter <- 0L
for (sp in seq_along(species_keys)) {
  key <- species_keys[sp]
  path <- file.path(input_dir,paste0(key,".regional_plot_data.tsv.gz"))
  if (!file.exists(path)) stop("Missing original input file: ",path)
  con <- gzfile(path,"rt")
  dat <- tryCatch(read.delim(con,stringsAsFactors=FALSE,check.names=FALSE),finally=close(con))
  stopifnot(nrow(dat)==expected_n[sp])
  for (strand in c("plus","minus")) {
    cpm <- dat[[paste0("mean_cpm_",strand)]]
    stopifnot(all(is.finite(cpm)),all(cpm>=0))
    x <- log2(1+cpm)
    stopifnot(max(x)<X_MAX)
    # Exact left-closed bins, matching numpy.histogram. Zero belongs to the first bin.
    bin_id <- findInterval(x,edges,rightmost.closed=TRUE)
    stopifnot(all(bin_id>=1L),all(bin_id<length(edges)))
    counts <- tabulate(bin_id,nbins=length(edges)-1L)
    frequency <- 100*counts/length(x)
    stopifnot(sum(counts)==length(x),abs(sum(frequency)-100)<1e-10)
    panel <- list(counts=counts,frequency=frequency,n=length(x),zero=sum(cpm==0),
                  zero_percent=100*mean(cpm==0))
    panels[[paste(key,strand,sep=".")]] <- panel
    for (j in seq_along(counts)) {
      counter <- counter+1L
      bin_rows[[counter]] <- data.frame(species=key,strand=strand,bin_index=j,
        left_log2p1_cpm=edges[j],right_log2p1_cpm=edges[j+1],window_count=counts[j],
        frequency_percent=frequency[j],n_total=length(x),stringsAsFactors=FALSE)
    }
  }
}

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
frequency_table <- do.call(rbind,bin_rows)
reference <- read.csv(file.path(figure_dir,"data/frequency/true_signal_frequency_bins.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
stopifnot(nrow(frequency_table)==nrow(reference))
for(column in names(frequency_table)) {
  if(is.numeric(frequency_table[[column]])) {
    stopifnot(isTRUE(all.equal(frequency_table[[column]],reference[[column]],tolerance=1e-12)))
  } else stopifnot(identical(frequency_table[[column]],reference[[column]]))
}
write.csv(frequency_table,file.path(output_dir,"frequency_bins_R.csv"),row.names=FALSE)
mm <- function(x) grid::unit(x,"mm")
text_mm <- function(label,x,y,just="left",bold=FALSE,italic=FALSE,size=FONT_PT,color="#222222",rot=0) {
  face <- if (bold) "bold" else if (italic) "italic" else "plain"
  grid::grid.text(label,x=mm(x),y=mm(y),just=just,rot=rot,
    gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=size,fontface=face,col=color))
}
line_mm <- function(x,y,color="#222222",width=.5) {
  grid::grid.lines(x=mm(x),y=mm(y),gp=grid::gpar(col=color,lwd=width,lineend="butt"))
}
draw_panel <- function(sp,si) {
  key <- species_keys[sp];strand <- c("plus","minus")[si]
  p <- panels[[paste(key,strand,sep=".")]]
  left <- 13.5+15+(sp-1)*57.5
  bottom <- 106+if (si==1L) 89 else 15
  width <- 47.5;height <- 46
  y_position <- function(value) {
    if(Y_LOG) bottom+(log10(value)-log10(LOGY_LIMITS[1]))/diff(log10(LOGY_LIMITS))*height
    else bottom+value/100*height
  }
  color <- if (si==1L) PLUS_COLOR else MINUS_COLOR
  header <- 106+if (si==1L) 144 else 70
  text_mm(letters[(si-1L)*3+sp],left-6,header,bold=TRUE,size=PANEL_PT)
  text_mm(if (si==1L) "Plus strand (+)" else "Minus strand (-)",
          left+width/2,header,just="center",color=color)
  for (j in seq_along(p$frequency)) {
    if (p$frequency[j]>0) {
      grid::grid.rect(x=mm(left+edges[j]/X_MAX*width),y=mm(bottom),
        width=mm(BIN_WIDTH/X_MAX*width),height=mm(y_position(p$frequency[j])-bottom),
        just=c("left","bottom"),gp=grid::gpar(fill=color,col=NA,alpha=.9))
    }
  }
  line_mm(c(left,left+width),c(bottom,bottom))
  line_mm(c(left,left),c(bottom,bottom+height))
  for (v in if(Y_LOG) c(0.001,0.01,0.1,1,10,100) else seq(0,100,by=25)) {
    yy <- y_position(v)
    line_mm(c(left-.7,left),c(yy,yy))
    text_mm(as.character(v),left-1.2,yy,just="right")
  }
  for (v in c(0,3,6,9,12)) {
    xx <- left+v/X_MAX*width
    line_mm(c(xx,xx),c(bottom-.7,bottom))
    text_mm(as.character(v),xx,bottom-2.6,
            just=if (v==0) "left" else if (v==X_MAX) "right" else "center")
  }
  if (sp==1L) text_mm(if(Y_LOG) "Windows (%; log scale)" else "Windows (%)",left-7.5,bottom+height/2,just="center",rot=90)
  text_mm(paste0("N = ",format(p$n,big.mark=",",trim=TRUE)),left+width-1,bottom+height-2.7,just="right")
  text_mm(sprintf("CPM = 0: %s (%.1f%%)",format(p$zero,big.mark=",",trim=TRUE),p$zero_percent),
          left+width-1,bottom+height-7.2,just="right",color="#555555")
}
draw_page <- function(crop=FALSE) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x=mm(if(crop) -13.5 else 0),y=mm(if(crop) -108 else 0),width=mm(210),height=mm(297),just=c("left","bottom"),clip="off"))
  on.exit(grid::popViewport(),add=TRUE)
  text_mm(if(Y_LOG) "Distribution of experimental GRO-seq signal (logarithmic y axis)" else "Distribution of experimental GRO-seq signal",28.5,271,bold=TRUE)
  text_mm("strand_seiplant_v1 | all 1,024-bp test windows | repeats retained",28.5,266)
  for (sp in 1:3) {
    text_mm(scientific_names[sp],13.5+15+(sp-1)*57.5+23.75,257,just="center",italic=TRUE)
    draw_panel(sp,1L);draw_panel(sp,2L)
  }
  text_mm("Experimental signal, log2(1 + mean SRX CPM)",105,111,just="center")
  text_mm("Experimental GRO-seq signal frequency in held-out test windows.",13.5,96,bold=TRUE)
  caption <- c(
    "a-c, Plus strand; d-f, minus strand. Columns show Arabidopsis, wheat and maize. Histograms include all formal",
    "held-out test windows, including CPM = 0. The x axis is log2(1 + mean SRX CPM), using common 0.25-wide",
    "bins from 0 to 12. The y axis is each bin's count divided by all windows in that panel, multiplied by 100%.",
    "Frequencies sum to 100%. The first bin includes zero signals and low positive signals; the separate CPM = 0",
    "annotation gives the exact zero count and percentage. Each SRX was independently CPM-normalized before",
    "equal-weight averaging. Arabidopsis / wheat / maize: 45,066 / 98,656 / 49,786 windows; 14 / 5 / 3 SRX.",
    "The original five-library wheat results are used. Repeat sequences are retained. No windows were subsampled",
    "and no signal values were jittered or removed. Bins are left-closed and right-open, except the final bin.",
    "The distribution describes the supplied test-window sampling design rather than random whole-genome windows.",
    if(Y_LOG) "Logarithmic y detail: percentages and denominators are unchanged; empty bins are not drawn on the log axis."
      else "This R script independently recomputes frequencies from the supplied original window-level data."
  )
  for (i in seq_along(caption)) text_mm(caption[i],13.5,91.5-(i-1)*3.5)
}

for (Y_LOG in if(MAKE_LOG_Y) c(FALSE,TRUE) else FALSE) {
  if(Y_LOG) {
    positive <- frequency_table$frequency_percent[frequency_table$frequency_percent>0]
    stopifnot(min(positive)>=LOGY_LIMITS[1],max(positive)<=LOGY_LIMITS[2])
  }

for (crop in c(FALSE, TRUE)) {
  dimensions <- if(crop) c(183,165) else c(210,297)
  stem <- paste0("true_signal_frequency",if(Y_LOG) "_logY_detail" else "",if(crop) "" else "_A4")
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
}
