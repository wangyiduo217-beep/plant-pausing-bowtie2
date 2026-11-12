# Run this file with Source in RStudio; edit config.R to compare visual styles.
# No CRAN packages, Python, network, or external data paths are needed.
find_bundle_dir <- function() {
  override <- Sys.getenv("GROSEQ_BUNDLE_DIR", "")
  if (nzchar(override)) return(normalizePath(override, winslash="/", mustWork=TRUE))
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f) && length(f) == 1L && file.exists(f)) {
      d <- dirname(normalizePath(f, winslash="/", mustWork=TRUE))
      if (file.exists(file.path(d,"config.R")) && file.exists(file.path(d,"run_plots.R"))) return(d)
    }
  }
  args <- commandArgs(trailingOnly=FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE)))
  if (file.exists(file.path(getwd(),"config.R")) && file.exists(file.path(getwd(),"run_plots.R"))) return(getwd())
  stop("Cannot locate data/. Open the .Rproj file, or use source('full/path/run_plots.R').")
}
bundle_dir <- find_bundle_dir()
figure_dir <- normalizePath(file.path(bundle_dir,"..",".."),winslash="/",mustWork=TRUE)
data_dir <- file.path(figure_dir,"data","regional")
track_data_dir <- file.path(figure_dir,"data","tracks")
config_file <- Sys.getenv("GROSEQ_CONFIG",file.path(bundle_dir,"config.R"))
sys.source(config_file, envir=environment())
threshold_override <- Sys.getenv("GROSEQ_MIN_PREDICTION_SCORE", "")
if (nzchar(threshold_override)) MIN_PREDICTION_SCORE <- as.numeric(threshold_override)
stopifnot(POINT_STYLE %in% c("density","scatter","hexbin"), FIT_METHOD %in% c("ols","local_linear","none"),
          Y_TRANSFORM %in% c("log2p1","log2","raw"), is.numeric(X_LIMITS), length(X_LIMITS)==2,
          length(MIN_PREDICTION_SCORE)==1L, is.finite(MIN_PREDICTION_SCORE),
          MIN_PREDICTION_SCORE>=0, MIN_PREDICTION_SCORE<=1)
if (!capabilities("cairo")) stop("This R build needs Cairo graphics support (standard Windows R supports it).")
threshold_tag <- format(MIN_PREDICTION_SCORE, scientific=FALSE, trim=TRUE, digits=12)
run_tag <- paste(POINT_STYLE,FIT_METHOD,Y_TRANSFORM,paste0("min_score_",threshold_tag),sep="_")
output_dir <- file.path(figure_dir,"results",OUTPUT_SUBDIR,run_tag)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
spec <- data.frame(key=c("arabidopsis","wheat","maize"),
    scientific=c("Arabidopsis thaliana","Triticum aestivum","Zea mays"),
    chrom=c("3","3A / 3B / 3D","6"), n=c(45066L,98656L,49786L),srx=c(14L,5L,3L),
    stringsAsFactors=FALSE)
transform_y <- function(score) {
  if (Y_TRANSFORM=="log2p1") return(log2(1+score))
  if (Y_TRANSFORM=="log2") { if(any(score<=0)) stop("log2 requires all scores > 0.");return(log2(score)) }
  score
}
y_label <- switch(Y_TRANSFORM,log2p1="log2(1 + predicted score)",log2="Predicted score (log2)",raw="Predicted score")
if (is.null(Y_LIMITS)) Y_LIMITS <- if(Y_TRANSFORM=="log2") c(-20.2,0.2) else c(-0.1,1.02)
y_ticks <- if(Y_TRANSFORM=="log2") c(-20,-10,0) else c(0,0.5,1)
x_ticks <- c(0,6,12)
x_label <- "GRO-seq signal (log2[1 + mean CPM])"

# A binned Gaussian density estimate changes colors only, never point positions.
# All operations below use base R; the matrix is 256 x 256, not an n x n KDE.
relative_density <- function(x,y) {
  nx <- 256L; ny <- 256L; xr <- c(0,12); yr <- if(Y_TRANSFORM=="log2") c(-20,0) else c(0,1)
  ix <- pmax(1L,pmin(nx,floor((x-xr[1])/diff(xr)*nx)+1L))
  iy <- pmax(1L,pmin(ny,floor((y-yr[1])/diff(yr)*ny)+1L))
  h <- matrix(tabulate((ix-1L)*ny+iy,nbins=nx*ny),nrow=ny,ncol=nx)
  stopifnot(sum(h)==length(x))
  reflect <- function(index,n) {index[index<1L] <- 1L-index[index<1L];index[index>n] <- 2L*n+1L-index[index>n];index}
  offsets <- -8:8; weights <- exp(-0.5*(offsets/2)^2);weights <- weights/sum(weights)
  hy <- matrix(0,ny,nx)
  for(k in seq_along(offsets)) hy <- hy+weights[k]*h[reflect(seq_len(ny)+offsets[k],ny),,drop=FALSE]
  hs <- matrix(0,ny,nx)
  for(k in seq_along(offsets)) hs <- hs+weights[k]*hy[,reflect(seq_len(nx)+offsets[k],nx),drop=FALSE]
  tx <- pmax(0,pmin(nx-1,(x-xr[1])/diff(xr)*nx-0.5))
  ty <- pmax(0,pmin(ny-1,(y-yr[1])/diff(yr)*ny-0.5))
  x0 <- floor(tx)+1L;y0 <- floor(ty)+1L;x1 <- pmin(nx,x0+1L);y1 <- pmin(ny,y0+1L)
  fx <- tx-floor(tx);fy <- ty-floor(ty)
  d <- (1-fx)*(1-fy)*hs[cbind(y0,x0)]+fx*(1-fy)*hs[cbind(y0,x1)]+
        (1-fx)*fy*hs[cbind(y1,x0)]+fx*fy*hs[cbind(y1,x1)]
  log1p(d)/log1p(max(d))
}
make_hex_data <- function(x,y) {
  # Two interleaved lattices, with nearest-center allocation in the hexagon metric.
  # Correlations and OLS are still calculated from all original observations.
  nx <- as.integer(HEX_GRIDSIZE[1]);ny <- as.integer(HEX_GRIDSIZE[2])
  xr <- c(-12e-9,12+12e-9);yr <- if(Y_TRANSFORM=="log2")c(-20,0) else c(0,1)
  sx <- diff(xr)/nx;sy <- diff(yr)/ny
  tx <- (x-xr[1])/sx;ty <- (y-yr[1])/sy
  ix1 <- round(tx);iy1 <- round(ty);ix2 <- floor(tx);iy2 <- floor(ty)
  use1 <- (tx-ix1)^2+3*(ty-iy1)^2 < (tx-ix2-0.5)^2+3*(ty-iy2-0.5)^2
  n1 <- (nx+1L)*(ny+1L)
  ids <- ifelse(use1,ix1*(ny+1L)+iy1+1L,n1+ix2*ny+iy2+1L)
  stopifnot(all(ids>=1L & ids<=n1+nx*ny))
  counts <- tabulate(ids,nbins=n1+nx*ny)
  stopifnot(sum(counts)==length(x))
  grid1 <- expand.grid(y=0:ny,x=0:nx)
  grid2 <- expand.grid(y=seq_len(ny)-0.5,x=seq_len(nx)-0.5)
  centers <- rbind(grid1,grid2);keep <- counts>0L
  list(x=xr[1]+centers$x[keep]*sx,y=yr[1]+centers$y[keep]*sy,
       counts=counts[keep],sx=sx,sy=sy,n_counted=sum(counts))
}

make_fit <- function(x,y) {
  xx <- seq(min(x),max(x),length.out=401)
  if(FIT_METHOD=="none") return(data.frame(x=xx,y=NA_real_,ess=NA_real_))
  if(FIT_METHOD=="ols") {
    model <- lm(y~x)
    return(data.frame(x=xx,y=unname(coef(model)[1]+coef(model)[2]*xx),ess=NA_real_))
  }
  yy <- ess <- numeric(length(xx))
  for(i in seq_along(xx)) {
    d <- x-xx[i];w <- exp(-0.5*(d/LOCAL_BANDWIDTH)^2)
    s0 <- sum(w);s1 <- sum(w*d);s2 <- sum(w*d*d);t0 <- sum(w*y);t1 <- sum(w*d*y)
    yy[i] <- (s2*t0-s1*t1)/(s0*s2-s1*s1);ess[i] <- s0*s0/sum(w*w)
  }
  yy[ess<MIN_EFFECTIVE_WEIGHT_COUNT] <- NA_real_
  data.frame(x=xx,y=yy,ess=ess)
}
panels <- list();metrics <- list();tracks <- list()
expected <- read.csv(file.path(data_dir,"expected_statistics_log2p1.csv"),stringsAsFactors=FALSE)
for (si in seq_len(nrow(spec))) {
  key <- spec$key[si]
  d <- read.delim(gzfile(file.path(data_dir,paste0(key,".regional_plot_data.tsv.gz"))),stringsAsFactors=FALSE)
  stopifnot(nrow(d)==spec$n[si],identical(d$window_index,seq_len(nrow(d))))
  for(strand in c("plus","minus")) {
    x <- log2(1+d[[paste0("mean_cpm_",strand)]])
    score <- d[[paste0("prediction_",strand)]];y <- transform_y(score)
    if(any(!is.finite(x)) || any(!is.finite(y))) stop("Nonfinite point detected; no records have been dropped.")
    full_n <- length(score)
    keep <- score >= MIN_PREDICTION_SCORE
    if(sum(keep)<3L) stop("Fewer than 3 retained windows for ", key, " / ", strand, "; lower MIN_PREDICTION_SCORE.")
    # The raw regression score controls selection, not a transformed or rounded value.
    original_index <- d$window_index[keep]
    removed_xpositive_n <- sum(!keep & x>0)
    removed_CPM_ge1_n <- sum(!keep & d[[paste0("mean_cpm_",strand)]]>=1)
    x <- x[keep];y <- y[keep];score <- score[keep]
    if(sd(x)==0 || sd(y)==0) stop("Retained data have zero variance for ", key, " / ", strand, "; correlation is undefined.")
    if(any(x<X_LIMITS[1] | x>X_LIMITS[2] | y<Y_LIMITS[1] | y>Y_LIMITS[2]))
      stop("Some points fall outside axis limits. Expand X_LIMITS/Y_LIMITS to retain all points.")
    r <- cor(x,y,method="pearson");rho <- cor(x,y,method="spearman");ols <- coef(lm(y~x))
    if(Y_TRANSFORM=="log2p1" && MIN_PREDICTION_SCORE==0) {
      ref <- expected[expected$species==key & expected$strand==strand,]
      stopifnot(abs(r-ref$pearson_r)<1e-10,abs(rho-ref$spearman_rho)<1e-10,
                abs(ols[1]-ref$OLS_intercept)<1e-10,abs(ols[2]-ref$OLS_slope)<1e-10)
    }
    density <- if(POINT_STYLE=="density") relative_density(x,y) else rep(0,length(x))
    colors <- if(POINT_STYLE=="density") colorRampPalette(DENSITY_PALETTE)(256)[1L+pmin(255L,floor(255*density))] else rep(SCATTER_COLOR,length(x))
    colors <- adjustcolor(colors,alpha.f=POINT_ALPHA)
    order <- if(POINT_STYLE=="density") order(density,seq_along(density)) else seq_along(x)
    fit <- make_fit(x,y);token <- paste(key,strand,sep=".")
    hex <- if(POINT_STYLE=="hexbin")make_hex_data(x,y) else NULL
    panels[[token]] <- list(x=x,y=y,colors=colors,order=order,fit=fit,r=r,rho=rho,n=length(x),
        original_n=full_n,original_index=original_index,hex=hex)
    metrics[[token]] <- data.frame(species=key,strand=strand,n=length(x),pearson_r=r,spearman_rho=rho,
        OLS_intercept=unname(ols[1]),OLS_slope=unname(ols[2]),transform_y=Y_TRANSFORM,fit_method=FIT_METHOD,
        threshold_raw_score=MIN_PREDICTION_SCORE,original_n=full_n,removed_n=full_n-length(x),
        retained_fraction=length(x)/full_n,removed_xpositive_n=removed_xpositive_n,
        removed_CPM_ge1_n=removed_CPM_ge1_n)
    write.csv(fit,file.path(output_dir,paste0(token,".fitted_line.csv")),row.names=FALSE)
    index_connection <- gzfile(file.path(output_dir,paste0(token,".retained_window_indices.csv.gz")),"wt")
    write.csv(data.frame(window_index=original_index),index_connection,row.names=FALSE)
    close(index_connection)
  }
  tracks[[key]] <- read.csv(file.path(track_data_dir,paste0(key,".selected_track.csv")),stringsAsFactors=FALSE)
}
metrics <- do.call(rbind,metrics)
write.csv(metrics,file.path(output_dir,"calculated_statistics.csv"),row.names=FALSE)

mm <- function(x) x/25.4
open_device <- function(file,width_mm,height_mm) {
  if(grepl("\\.pdf$",file)) cairo_pdf(file,width=mm(width_mm),height=mm(height_mm),family=FONT_FAMILY,pointsize=FONT_PT,onefile=TRUE)
  else png(file,width=mm(width_mm),height=mm(height_mm),units="in",res=PNG_DPI,type="cairo",family=FONT_FAMILY,pointsize=FONT_PT)
  par(ps=FONT_PT,cex=1,cex.axis=1,cex.lab=1,cex.main=1,las=1,family=FONT_FAMILY)
  plot.new()
}
page_context <- list(width=183,height=170,ox=0,oy=0)
label <- function(x,y,text,bold=FALSE,italic=FALSE,size=FONT_PT,just="left",color="#111111") {
  grid::grid.text(text,x=grid::unit(x+page_context$ox,"mm"),y=grid::unit(y+page_context$oy,"mm"),just=just,
    gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=size,fontface=if(bold) "bold" else if(italic) "italic" else "plain",col=color))
}
panel_window <- function(x,y,w,h,xlim,ylim) {
  par(fig=c((x+page_context$ox)/page_context$width,(x+w+page_context$ox)/page_context$width,
            (y+page_context$oy)/page_context$height,(y+h+page_context$oy)/page_context$height),
      mai=c(0,0,0,0),new=TRUE,ps=FONT_PT,cex=1,cex.axis=1,cex.lab=1,las=1,xaxs="i",yaxs="i",family=FONT_FAMILY)
  plot.new();plot.window(xlim=xlim,ylim=ylim,xaxs="i",yaxs="i")
}
draw_scatter <- function(p,x,y,w,h,ylab=FALSE) {
  panel_window(x,y,w,h,X_LIMITS,Y_LIMITS)
  if(POINT_STYLE=="hexbin") {
    hx <- p$hex
    cols <- colorRampPalette(DENSITY_PALETTE)(256)[1L+pmin(255L,floor(255*log1p(hx$counts)/log1p(max(hx$counts))))]
    vx <- c(0.5,0.5,0,-0.5,-0.5,0)*hx$sx
    vy <- c(-0.5,0.5,1,0.5,-0.5,-1)*hx$sy/3
    for(i in seq_along(hx$counts)) polygon(hx$x[i]+vx,hx$y[i]+vy,col=cols[i],border=NA)
  } else {
    o <- p$order;points(p$x[o],p$y[o],pch=16,cex=POINT_CEX,col=p$colors[o])
  }
  if(FIT_METHOD!="none") lines(p$fit$x,p$fit$y,col=FIT_COLOR,lwd=FIT_LINE_WIDTH,lty=FIT_LINE_TYPE)
  axis(1,at=x_ticks,labels=format(x_ticks,trim=TRUE),tck=-0.025,mgp=c(0,0.33,0),lwd=0.5,lwd.ticks=0.5)
  axis(2,at=y_ticks,labels=format(y_ticks,trim=TRUE),tck=-0.025,mgp=c(0,0.33,0),lwd=0.5,lwd.ticks=0.5,las=1)
  if(ylab) grid::grid.text(y_label,x=grid::unit(x-8+page_context$ox,"mm"),
      y=grid::unit(y+h/2+page_context$oy,"mm"),rot=90,
      gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT))
}
draw_bottom <- function(full=FALSE) {
  label(if(full)0 else 18,2.5,sprintf("strand_seiplant_v1; raw score >= %s; retained windows only; 14 / 5 / 3 SRX",threshold_tag))
  xx <- if(full)0 else 18
  grid::grid.lines(x=grid::unit(c(xx,xx+7)+page_context$ox,"mm"),y=grid::unit(c(9,9)+page_context$oy,"mm"),
      gp=grid::gpar(col=FIT_COLOR,lwd=FIT_LINE_WIDTH,lty=FIT_LINE_TYPE))
  label(xx+9,9,if(FIT_METHOD=="ols")"OLS fit" else if(FIT_METHOD=="local_linear")"Local linear fit" else "No fit")
  if(full && any(vapply(tracks,function(t)any(!is.finite(t$prediction_plus) | !is.finite(t$prediction_minus)),logical(1)))) {
    grid::grid.rect(x=grid::unit(54+page_context$ox,"mm"),y=grid::unit(9+page_context$oy,"mm"),
      width=grid::unit(3,"mm"),height=grid::unit(1.8,"mm"),gp=grid::gpar(fill="#e6e6e6",col=NA))
    label(57,9,"Missing prediction")
  }
  legend_name <- switch(POINT_STYLE,density="Within-panel density",hexbin="Relative hexagon count",scatter="All retained windows")
  label(119,9,legend_name,just="right")
  if(POINT_STYLE %in% c("density","hexbin")) {
    ramp <- colorRampPalette(DENSITY_PALETTE)(100)
    grid::grid.raster(as.raster(matrix(ramp,nrow=1)),x=grid::unit(153+page_context$ox,"mm"),
      y=grid::unit(9+page_context$oy,"mm"),width=grid::unit(57,"mm"),height=grid::unit(1.6,"mm"),interpolate=TRUE)
    label(124.5,5.8,"Low");label(181.5,5.8,"High",just="right")
  }
}
draw_six <- function() {
  label(18,167,sprintf("GRO-seq predictions | raw score >= %s",threshold_tag),bold=TRUE)
  for(si in seq_len(nrow(spec))) {
    x <- 18+58*(si-1);mid <- x+23.5
    label(mid,160,spec$scientific[si],italic=TRUE,just="center")
    label(mid,155,sprintf("chr%s",spec$chrom[si]),just="center")
    for(j in 1:2) {
      strand <- c("plus","minus")[j];p <- panels[[paste(spec$key[si],strand,sep=".")]];top <- 148-(j-1)*66
      label(x-7,top,letters[(j-1)*3+si],bold=TRUE,size=PANEL_PT)
      label(mid,top,sprintf("%s; n = %s",if(j==1)"Plus strand" else "Minus strand",format(p$n,big.mark=",",trim=TRUE)),just="center")
      label(mid,top-4,sprintf("r = %.3f; \u03c1 = %.3f",p$r,p$rho),just="center",color=STAT_COLOR)
      draw_scatter(p,x,top-53,47,46,ylab=(si==1))
    }
  }
  label(99.5,18,x_label,just="center");draw_bottom()
}
regions <- read.csv(file.path(track_data_dir,"selected_regions.csv"),stringsAsFactors=FALSE)
draw_composite <- function() {
  label(21,167,"Representative genomic region",bold=TRUE);label(113,167,sprintf("Test windows: raw score >= %s",threshold_tag),bold=TRUE)
  track_colors <- c("#21844B","#78B799","#3278A5","#8BAFD2")
  track_names <- c("GRO-seq +","Prediction +","GRO-seq -","Prediction -")
  for(si in seq_len(nrow(spec))) {
    key <- spec$key[si];t <- tracks[[key]];rr <- regions[regions$species==key,];top <- 160-(si-1)*51
    label(0,top,letters[(si-1)*3+1],bold=TRUE,size=PANEL_PT)
    label(21,top,spec$scientific[si],italic=TRUE)
    label(21,top-3.5,sprintf("chr%s:%s-%s",rr$chrom,format(rr$start0+1,big.mark=",",scientific=FALSE,trim=TRUE),format(rr$end0,big.mark=",",scientific=FALSE,trim=TRUE)))
    pn <- panels[[paste(key,"plus",sep=".")]]$n;mn <- panels[[paste(key,"minus",sep=".")]]$n
    label(147.75,top+3.5,sprintf("chr%s; kept +/- n = %s / %s",spec$chrom[si],format(pn,big.mark=",",trim=TRUE),format(mn,big.mark=",",trim=TRUE)),just="center")
    values <- list(t$gro_log2p1_plus,t$prediction_plus,t$gro_log2p1_minus,t$prediction_minus)
    for(k in 1:4) {
      low <- top-8-6.4-(k-1)*8.5;maxy <- if(k %in% c(2,4))1 else max(0.5,ceiling(max(values[[k]])*1.08*2)/2)
      panel_window(21,low,77,6.4,c(rr$start0,rr$end0)/1e6,c(0,maxy))
      missing <- !is.finite(values[[k]])
      if(k %in% c(2,4) && "valid" %in% names(t)) missing <- missing | !as.logical(t$valid)
      if(k %in% c(2,4) && any(missing)) rect(t$start0[missing]/1e6,0,t$end0[missing]/1e6,maxy,col="#e6e6e6",border=NA)
      keep <- !missing
      rect(t$start0[keep]/1e6,0,t$end0[keep]/1e6,values[[k]][keep],col=track_colors[k],border=NA)
      axis(2,at=c(0,maxy),labels=c("0",format(maxy,trim=TRUE)),tck=-0.02,mgp=c(0,0.33,0),las=1,lwd=0.5,lwd.ticks=0.5)
      tick_step <- if((rr$end0-rr$start0)>60000)20000 else 10000
      xticks <- seq(rr$start0+tick_step/2,rr$end0-1,by=tick_step)/1e6
      axis(1,at=xticks,labels=if(k==4)format(xticks,nsmall=3,trim=TRUE) else FALSE,
           tck=if(k==4)-0.12 else 0,mgp=c(0,0.33,0),lwd=0.5,lwd.ticks=0.5)
      label(14.3,low+3.2,track_names[k],just="right")
      if(k==4)label(59.5,low-5.2,sprintf("Position on chr%s (Mb)",rr$chrom),just="center")
    }
    for(j in 1:2) {
      strand <- c("plus","minus")[j];p <- panels[[paste(key,strand,sep=".")]];x <- 113+(j-1)*39
      label(x-9.5,top,letters[(si-1)*3+j+1],bold=TRUE,size=PANEL_PT)
      label(x+15.25,top,if(j==1)"Plus strand" else "Minus strand",just="center")
      label(x+15.25,top-3.5,sprintf("r = %.3f; \u03c1 = %.3f",p$r,p$rho),just="center",color=STAT_COLOR)
      draw_scatter(p,x,top-38,30.5,30,ylab=(j==1))
    }
  }
  label(147.75,12.8,x_label,just="center");draw_bottom(TRUE)
}
caption <- paste("A prespecified raw-prediction-score filter is applied separately to each strand of the formal held-out test windows.",
  paste0("Keep original score >= ",threshold_tag,"; threshold zero retains every input window."),
  "Original Arabidopsis / wheat / maize test sets contain 45,066 / 98,656 / 49,786 windows, with 14 / 5 / 3 SRX; repeats retained.",
  "Each panel reports its retained n; no additional sampling, x=0 filtering, rounding, or coordinate shifting is applied.",
  "The x-axis is log2(1 + mean SRX CPM); each SRX was independently normalized before equal-weight averaging.",
  paste0("The y-axis transformation is ",Y_TRANSFORM," (log2p1 means log2[1 + original predicted score])."),
  "r and rho are recomputed Pearson and Spearman correlations on the retained original coordinates.",
  paste0("The line method is ",FIT_METHOD,"; OLS uses all displayed points, not y=x."),
  paste0("Display style: ",POINT_STYLE,"; hexbin aggregates all windows while statistics still use the original observations."),
  "Both axes start at zero by default. Thresholds refer to regression scores, not probabilities or confidence levels.",
  "This is a conditional-subset display sensitivity comparison, not an improvement in the frozen model's full-test performance.",
  "In the composite, selected track examples retain the original GRO-seq log2(1+CPM) and 0-1 prediction scales. The shared track tables give each exact region.",
  "Light-gray prediction-track bins mark unavailable predictions (reference N input), not zero scores; the corresponding experimental signal remains shown.")
writeLines(caption,file.path(output_dir,"figure_legend_EN.txt"),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(output_dir,"R_sessionInfo.txt"))
writeLines(readLines(config_file,warn=FALSE),file.path(output_dir,"settings_used.R"),useBytes=TRUE)
write.csv(data.frame(axis=c("x","y"),min=c(X_LIMITS[1],Y_LIMITS[1]),max=c(X_LIMITS[2],Y_LIMITS[2])),
          file.path(output_dir,"axis_limits_used.csv"),row.names=FALSE)
export <- function(stem,draw,width=183,height=170,a4=FALSE) {
  for(ext in c("pdf","png")) {
    file <- file.path(output_dir,paste0(stem,".",ext));open_device(file,if(a4)210 else width,if(a4)297 else height)
    page_context <<- list(width=if(a4)210 else width,height=if(a4)297 else height,ox=if(a4)13.5 else 0,oy=if(a4)106 else 0)
    draw()
    if(a4) {
      grid::grid.text("GRO-seq signal and sequence-based predictions",x=grid::unit(13.5,"mm"),y=grid::unit(96,"mm"),just="left",gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT,fontface="bold"))
      for(i in seq_along(strwrap(caption,width=152)))
        grid::grid.text(strwrap(caption,width=152)[i],x=grid::unit(13.5,"mm"),y=grid::unit(92-(i-1)*3.5,"mm"),just="left",gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT))
    }
    dev.off();message("Created: ",file)
  }
}
if(DRAW_CORRELATIONS) {
  export("three_species_correlations_183x170mm",draw_six)
  if(MAKE_A4_PROOFS)export("three_species_correlations_A4",draw_six,a4=TRUE)
}
if(MAKE_TRACK_COMPOSITE) {
  export("tracks_and_correlations_183x170mm",draw_composite)
  if(MAKE_A4_PROOFS)export("tracks_and_correlations_A4",draw_composite,a4=TRUE)
}
if(MAKE_INDIVIDUAL_SPECIES) for(si in seq_len(nrow(spec))) {
  local({idx <- si
    draw_pair <- function() {
      label(51,79,spec$scientific[idx],italic=TRUE,just="center")
      for(j in 1:2) {
        x <- 12+(j-1)*47;p <- panels[[paste(spec$key[idx],c("plus","minus")[j],sep=".")]]
        label(x-6,71,letters[j],bold=TRUE,size=PANEL_PT)
        label(x+17,71,sprintf("%s; n = %s",if(j==1)"Plus strand" else "Minus strand",format(p$n,big.mark=",",trim=TRUE)),just="center")
        label(x+17,67,sprintf("r = %.3f; \u03c1 = %.3f",p$r,p$rho),just="center",color=STAT_COLOR)
        draw_scatter(p,x,20,34,43,ylab=(j==1))
      }
      label(52,8,x_label,just="center")
    }
    export(paste0(spec$key[idx],"_two_strands"),draw_pair,width=100,height=85)
  })
}
message("Finished. Results are in: ",output_dir)
