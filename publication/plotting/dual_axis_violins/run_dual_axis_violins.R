# Standard R only. Shared data are in publication/data/regional, located relative to this script.
find_dir <- function() {
  p <- Sys.getenv("GROSEQ_DUAL_VIOLIN_DIR", "")
  if(nzchar(p))return(normalizePath(p,winslash="/",mustWork=TRUE))
  for(fr in rev(sys.frames()))if(!is.null(fr$ofile)&&file.exists(fr$ofile)) {
    p <- dirname(normalizePath(fr$ofile,winslash="/",mustWork=TRUE))
    if(file.exists(file.path(p,"config.R"))&&file.exists(file.path(p,"run_dual_axis_violins.R")))return(p)
  }
  f <- commandArgs(FALSE);f <- sub("^--file=","",f[grepl("^--file=",f)])
  if(length(f))return(dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE)))
  if(file.exists("run_dual_axis_violins.R")&&file.exists("config.R"))return(getwd())
  stop("Source the full path to run_dual_axis_violins.R or open the supplied project.")
}
bundle <- find_dir();sys.source(file.path(bundle,"config.R"),envir=environment())
stopifnot(capabilities("cairo"),EXPERIMENT_Y_MAX>0,PREDICTION_Y_MAX>0)
publication_dir <- normalizePath(file.path(bundle,"..",".."),winslash="/",mustWork=TRUE)
data_dir <- normalizePath(file.path(publication_dir,"data","regional"),winslash="/",mustWork=TRUE)
out <- file.path(publication_dir,"generated","dual_axis_violins");dir.create(out,recursive=TRUE,showWarnings=FALSE)
keys <- c("arabidopsis","wheat","maize")
scientific <- c("Arabidopsis thaliana","Triticum aestivum","Zea mays")
regions <- c("upstream_1_3kb","upstream_lt1kb","five_prime_utr","CDS","intron","three_prime_utr","downstream_le1kb","distal")
labs <- c("Upstream\n(1-3 kb)","Upstream\n(<1 kb)","5' UTR","CDS","Intron","3' UTR","Downstream\n(<=1 kb)","Distal")
expected_n <- c(45066L,98656L,49786L);chrom <- c("3","3A / 3B / 3D","6");srx <- c(14L,5L,3L)
summary_vector <- function(v) {
  q <- unname(quantile(v,c(0,.25,.5,.75,.95,.99,1),type=7))
  setNames(c(q,mean(v),sd(v),q[4]-q[2]),c("min","q25","median","q75","p95","p99","max","mean","sd","IQR"))
}
# Every observation enters the full-range histogram. The display upper limit is
# never used here. Reject an undersized KDE range instead of clamping high values.
violin_shape <- function(v,upper) {
  stopifnot(length(v)>1L,all(is.finite(v)),all(v>=0),all(v<=upper))
  nb <- 512L;step <- upper/nb
  bin <- pmin(nb,pmax(1L,floor(v/step)+1L));h <- tabulate(bin,nbins=nb)
  stopifnot(sum(h)==length(v))
  iq <- IQR(v,type=7);spread <- if(iq>0)min(sd(v),iq/1.34)else sd(v)
  bw <- max(.9*spread*length(v)^(-.2)*BANDWIDTH_ADJUST,step*.7)
  sigma <- bw/step;radius <- as.integer(4*sigma+.5);delta <- seq.int(-radius,radius)
  weights <- exp(-.5*(delta/sigma)^2);weights <- weights/sum(weights)
  reflect <- function(k) {
    while(any(k<1L|k>nb)){k[k<1L] <- 1L-k[k<1L];k[k>nb] <- 2L*nb+1L-k[k>nb]}
    k
  }
  hs <- numeric(nb)
  for(k in seq_along(delta))hs <- hs+weights[k]*h[reflect(seq_len(nb)+delta[k])]
  yy <- seq(min(v),max(v),length.out=256L)
  den <- approx((seq_len(nb)-.5)*step,hs,xout=yy,rule=2,ties="ordered")$y
  list(y=yy,width=den/max(den)*VIOLIN_HALF_WIDTH,q=unname(quantile(v,c(.25,.5,.75),type=7)),
       n=length(v),bandwidth=bw,bin_min=0,bin_max=upper,bin_count=nb,histogram_count=sum(h))
}
datalist <- list();summaries <- list();kde_records <- list();shape_records <- list()
for(si in seq_along(keys)) {
  key <- keys[si]
  d <- read.delim(gzfile(file.path(data_dir,paste0(key,".regional_plot_data.tsv.gz"))),stringsAsFactors=FALSE)
  stopifnot(nrow(d)==expected_n[si],all(d$region%in%regions),all(d$window_index==seq_len(nrow(d))))
  dd <- list(region=d$region,counts=tabulate(match(d$region,regions),nbins=8L))
  for(strand in c("plus","minus"))for(kind in c("experimental","prediction")) {
    raw <- d[[paste0(if(kind=="experimental")"mean_cpm_"else"prediction_",strand)]]
    v <- log2(1+raw);stopifnot(all(is.finite(v)),all(raw>=0))
    if(kind=="prediction")stopifnot(all(raw<=1))
    upper <- if(kind=="experimental")EXPERIMENT_KDE_MAX else PREDICTION_KDE_MAX
    display <- if(kind=="experimental")EXPERIMENT_Y_MAX else PREDICTION_Y_MAX
    shapes <- list()
    for(j in seq_along(regions)) {
      take <- d$region==regions[j];rr <- raw[take];vv <- v[take]
      a <- violin_shape(vv,upper);shapes[[regions[j]]] <- a
      z <- data.frame(species=key,species_name=scientific[si],strand=strand,quantity=kind,region=regions[j],
        n_windows=length(vv),raw_zero_count=sum(rr==0),raw_zero_fraction=mean(rr==0),display_y_max=display,
        above_display_count=sum(vv>display),above_display_fraction=mean(vv>display),stringsAsFactors=FALSE)
      for(prefix in c("raw","log2p1")) {
        st <- summary_vector(if(prefix=="raw")rr else vv)
        for(nm in names(st))z[[paste0(prefix,"_",nm)]] <- unname(st[[nm]])
      }
      summaries[[length(summaries)+1L]] <- z
      kde_records[[length(kde_records)+1L]] <- data.frame(species=key,strand=strand,quantity=kind,region=regions[j],
        n_input=a$n,histogram_count=a$histogram_count,kde_bin_min=0,kde_bin_max=upper,kde_bin_count=a$bin_count,
        bandwidth=a$bandwidth,observed_min=min(vv),observed_max=max(vv),display_y_max=display,
        above_display_count=sum(vv>display),stringsAsFactors=FALSE)
      shape_records[[length(shape_records)+1L]] <- data.frame(species=key,strand=strand,quantity=kind,region=regions[j],
        support_index=seq_along(a$y),y_log2p1=a$y,half_width=a$width,stringsAsFactors=FALSE)
    }
    dd[[paste(strand,kind,sep=".")]] <- list(raw=raw,v=v,shapes=shapes)
  }
  datalist[[key]] <- dd
}
summary_table <- do.call(rbind,summaries);kde_table <- do.call(rbind,kde_records)
shape_table <- do.call(rbind,shape_records)
stopifnot(nrow(summary_table)==96L,nrow(kde_table)==96L,sum(expected_n)==193508L,all(kde_table$n_input==kde_table$histogram_count))
options(digits=17)
write.table(summary_table,file.path(out,"recalculated_signal_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(kde_table,file.path(out,"full_support_kde_parameters.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(shape_table,file.path(out,"full_support_violin_shapes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

ctx <- list(width=183,height=170,ox=0,oy=0);PX <- 22.5;PW <- 137;XR <- c(-.6,7.6)
txt <- function(x,y,label,size=FONT_PT,bold=FALSE,italic=FALSE,ha="left",col="#111111") {
  grid::pushViewport(grid::viewport(clip="off"))
  grid::grid.text(label,x=grid::unit(x+ctx$ox,"mm"),y=grid::unit(y+ctx$oy,"mm"),just=ha,
    gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=size,fontface=if(bold)"bold"else if(italic)"italic"else"plain",col=col))
  grid::popViewport()
}
rotated_label <- function(x,y,label,col) {
  grid::pushViewport(grid::viewport(clip="off"))
  grid::grid.text(label,x=grid::unit(x+ctx$ox,"mm"),y=grid::unit(y+ctx$oy,"mm"),rot=90,
    gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT,col=col));grid::popViewport()
}
group_x <- function(j)PX+(j-1-XR[1])/diff(XR)*PW
draw_violin <- function(a,pos,col) {
  # Clip only the final geometry; a contains the full observed support.
  par(xpd=FALSE)
  polygon(c(pos-a$width,rev(pos+a$width)),c(a$y,rev(a$y)),
    col=adjustcolor(col,.30),border=adjustcolor(col,.75),lwd=.55)
  segments(pos,a$q[1],pos,a$q[3],col=col,lwd=1.6)
  points(pos,a$q[2],pch=21,bg="white",col=col,cex=.45,lwd=.5,xpd=FALSE)
}
draw_panel <- function(si,strand,y,h=44) {
  d <- datalist[[keys[si]]]
  par(fig=c((PX+ctx$ox)/ctx$width,(PX+PW+ctx$ox)/ctx$width,(y+ctx$oy)/ctx$height,(y+h+ctx$oy)/ctx$height),
    mai=c(0,0,0,0),new=TRUE,ps=FONT_PT,cex=1,cex.axis=1,cex.lab=1,family=FONT_FAMILY,xaxs="i",yaxs="i",xpd=FALSE)
  plot.new();plot.window(XR,c(0,EXPERIMENT_Y_MAX),xaxs="i",yaxs="i");left_usr <- par("usr")
  axis(2,at=seq(0,EXPERIMENT_Y_MAX,length.out=5L),las=1,tck=-.015,lwd=.5,mgp=c(0,.3,0),cex.axis=1,
    col=EXPERIMENT_COLOR,col.axis=EXPERIMENT_COLOR)
  abline(h=seq(0,EXPERIMENT_Y_MAX,length.out=5L)[-c(1L,5L)],col="#EEEEEE",lwd=.35)
  segments(XR[1],0,XR[2],0,col="#888888",lwd=.5)
  for(j in seq_along(regions))draw_violin(d[[paste(strand,"experimental",sep=".")]]$shapes[[regions[j]]],j-1-.17,EXPERIMENT_COLOR)
  # Change native user coordinates for the right axis. Data stay in original
  # log2(1+score) units; there is NO min-max or per-region normalization.
  par(usr=c(XR,0,PREDICTION_Y_MAX),xpd=FALSE)
  stopifnot(isTRUE(all.equal(unname(par("usr")),c(XR,0,PREDICTION_Y_MAX))))
  axis(4,at=seq(0,PREDICTION_Y_MAX,length.out=5L),las=1,tck=-.015,lwd=.5,mgp=c(0,.3,0),cex.axis=1,
    col=PREDICTION_COLOR,col.axis=PREDICTION_COLOR)
  for(j in seq_along(regions))draw_violin(d[[paste(strand,"prediction",sep=".")]]$shapes[[regions[j]]],j-1+.17,PREDICTION_COLOR)
  par(usr=left_usr,xpd=FALSE)
  rotated_label(8.5,y+h/2,"GRO-seq: log2(1 + mean CPM)",EXPERIMENT_COLOR)
  rotated_label(175,y+h/2,"Prediction: log2(1 + score)",PREDICTION_COLOR)
  for(j in seq_along(regions)) {
    take <- d$region==regions[j]
    above <- mean(d[[paste(strand,"experimental",sep=".")]]$v[take]>EXPERIMENT_Y_MAX)*100
    txt(group_x(j),y+h+3.1,sprintf("\u2191 %.1f%%",above),ha="center",col=EXPERIMENT_COLOR)
    txt(group_x(j),y-3.5,paste0(labs[j],"\nn = ",format(d$counts[j],big.mark=",",trim=TRUE)),ha="top")
  }
}
legend_key <- function() {
  for(i in 1:2) {
    xx <- c(89,123)[i];col <- c(EXPERIMENT_COLOR,PREDICTION_COLOR)[i]
    grid::grid.rect(x=grid::unit(xx+ctx$ox+1.65,"mm"),y=grid::unit(159.5+ctx$oy,"mm"),
      width=grid::unit(3.3,"mm"),height=grid::unit(1.4,"mm"),gp=grid::gpar(fill=adjustcolor(col,.3),col=col,lwd=.5))
    txt(xx+4.7,159.5,c("Experiment","Prediction")[i])
  }
}
caption <- function(si) {
  title <- paste0(scientific[si],": regional experimental and predicted distributions.")
  lines <- c(
    "Green violins use the left axis, log2(1 + mean SRX CPM); orange violins use the right axis, log2(1 + predicted score).",
    "Two y axes have independent scales: compare distribution patterns, not absolute heights. No min-max normalization is used.",
    paste0("All ",format(expected_n[si],big.mark=",")," original held-out windows are retained, including zeros; no score threshold or outlier filtering is applied."),
    paste0("The left display ends at ",EXPERIMENT_Y_MAX,". Each upward-arrow percentage reports experimental values above this limit in that region."),
    "All values enter the full-support Gaussian density estimate before display clipping. Tails above the displayed limit are not deleted.",
    "White circles: medians; vertical bars: interquartile ranges. Equal maximum violin widths do not encode the number of windows.",
    "n is the full regional window count, not biological replicates. Regions use window midpoints and canonical protein-coding transcripts.",
    paste0("strand_seiplant_v1; repeats retained; ",srx[si]," SRX. Predictions are scores, not predicted CPM. No significance tests are performed."))
  grid::pushViewport(grid::viewport(clip="off"))
  grid::grid.text(title,x=grid::unit(13.5,"mm"),y=grid::unit(96,"mm"),just="left",gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT,fontface="bold"))
  for(i in seq_along(lines))grid::grid.text(lines[i],x=grid::unit(13.5,"mm"),y=grid::unit(91.9-(i-1)*3.5,"mm"),just="left",gp=grid::gpar(fontfamily=FONT_FAMILY,fontsize=FONT_PT))
  grid::popViewport()
}
draw_current <- function(si,a4=FALSE) {
  txt(PX,167,scientific[si],italic=TRUE)
  txt(159.5,167,"Regional GRO-seq and prediction distributions",ha="right",bold=TRUE)
  txt(PX,159.5,"Independent y axes; log2(1 + quantity)");legend_key()
  for(st in 1:2) {
    top <- c(150.5,77.5)[st]
    txt(1.5,top,letters[2*(si-1)+st],size=PANEL_PT,bold=TRUE)
    txt(PX,top,c("Plus strand (+)","Minus strand (-)")[st],bold=TRUE)
    txt(159.5,top,paste0("chr",chrom[si],"; n = ",format(expected_n[si],big.mark=","),"; ",srx[si]," SRX"),ha="right")
    draw_panel(si,c("plus","minus")[st],c(97,24)[st])
  }
  txt(PX,3.5,"Independent y axes: compare distribution patterns, not absolute heights. All test windows retained.")
  if(a4)caption(si)
}
open_page <- function(file,a4=FALSE) {
  pw <- if(a4)210 else 183;ph <- if(a4)297 else 170
  if(grepl("\\.pdf$",file))cairo_pdf(file,width=pw/25.4,height=ph/25.4,pointsize=FONT_PT,family=FONT_FAMILY,onefile=TRUE)
  else png(file,width=pw/25.4,height=ph/25.4,units="in",res=PNG_DPI,type="cairo",pointsize=FONT_PT,family=FONT_FAMILY)
  ctx <<- list(width=pw,height=ph,ox=if(a4)13.5 else 0,oy=if(a4)106 else 0)
  par(ps=FONT_PT,cex=1,cex.axis=1,cex.lab=1,family=FONT_FAMILY);plot.new()
}
for(a4 in c(FALSE,TRUE)) {
  if((a4&&!MAKE_A4)||(!a4&&!MAKE_MAIN))next
  suffix <- if(a4)"A4"else"183x170mm"
  file <- file.path(out,paste0("regional_dual_axis_",suffix,".pdf"));open_page(file,a4)
  for(si in 1:3){if(si>1){par(fig=c(0,1,0,1),new=FALSE);plot.new()};draw_current(si,a4)}
  dev.off();cat("Created:",file,"\n")
  for(si in 1:3) {
    file <- file.path(out,paste0(keys[si],"_dual_axis_",suffix,".png"))
    open_page(file,a4);draw_current(si,a4);dev.off()
  }
}
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
cat("Validated: all 193508 original windows; 96 distributions; full KDE; native dual axes; no filtering or min-max scaling.\n")
