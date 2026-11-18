# Portable, data-driven vector figure; all plotted scores retain their raw scale.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  for (loc in c("English_United States.utf8", "Chinese_China.utf8")) {
    ans <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error=function(e) ""))
    if (nzchar(ans)) break
  }
}
ROOT <- local({
  candidates <- character()
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) candidates <- c(candidates, dirname(normalizePath(f$ofile, winslash="/")))
  a <- commandArgs(FALSE); a <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(a)) candidates <- c(candidates, dirname(normalizePath(a[1], winslash="/")))
  candidates <- c(candidates, getwd())
  ok <- candidates[file.exists(file.path(candidates, "plot_radar_tracks_composite.R"))]
  if (!length(ok)) stop("Cannot locate the plotting script directory")
  ok[1]
})
library(grid)
# Repository wrapper only: keep the original numerical calculations and geometry.
args <- commandArgs(trailingOnly=TRUE)
output_arg <- which(args == "--output-dir")
if (length(output_arg) > 1L || (length(output_arg) && output_arg == length(args))) stop("Use --output-dir PATH once")
data_arg <- which(args == "--data-dir")
if (length(data_arg) > 1L || (length(data_arg) && data_arg == length(args))) stop("Use --data-dir PATH once")
standalone_layout <- dir.exists(file.path(ROOT,"data"))
DATA_DIR <- if (length(data_arg)) args[data_arg+1L] else if (standalone_layout) file.path(ROOT,"data") else file.path(ROOT,"..","..","data","radar_tracks")
DATA_DIR <- normalizePath(DATA_DIR,winslash="/",mustWork=TRUE)
DEFAULT_OUTPUT <- if (standalone_layout) file.path(ROOT,"output") else file.path(ROOT,"..","..","results","radar_tracks")
OUTPUT_DIR <- if (length(output_arg)) args[output_arg+1L] else DEFAULT_OUTPUT
DEFINE_ONLY <- "--define-only" %in% args || isTRUE(get0("GROSEQ_DEFINE_ONLY", ifnotfound=FALSE))
SOURCES <- c("arabidopsis", "maize", "wheat")
NAMES <- c(arabidopsis="Arabidopsis", maize="Maize", wheat="Wheat")
PCH <- c(arabidopsis=1, maize=2, wheat=0)
LTY <- c(arabidopsis="solid", maize="dashed", wheat="dotdash")
SAME <- "#FF7473"; CROSS <- "#34314C"
correlations <- read.csv(file.path(DATA_DIR,"correlation.csv"), stringsAsFactors=FALSE)
aurocs <- read.csv(file.path(DATA_DIR,"auroc.csv"), stringsAsFactors=FALSE)
aurocs <- subset(aurocs, strand == "micro")
peaks <- read.csv(file.path(DATA_DIR,"peak_overlap.csv"), stringsAsFactors=FALSE)
peaks <- subset(peaks, percentile == 85)
key <- c("target_species", "source_model")
stopifnot(nrow(correlations)==9, nrow(aurocs)==9, nrow(peaks)==9,
          all(correlations$strand=="pooled"), !anyDuplicated(correlations[key]),
          !anyDuplicated(aurocs[key]), !anyDuplicated(peaks[key]))
d <- merge(correlations[c(key, "n_observations", "pearson_r")],
           aurocs[c(key, "n_observations", "HiLearn_AUROC")], by=key, suffixes=c("_r", "_auc"))
d <- merge(d, peaks[c(key, "n_strand_window_pairs", "tp", "fp", "fn", "union_n", "micro_jaccard")], by=key)
stopifnot(nrow(d)==9, all(d$n_observations_r == d$n_observations_auc),
          all(d$n_observations_r == d$n_strand_window_pairs),
          all(abs(d$micro_jaccard - d$tp/d$union_n)<1e-12), all(d$union_n == d$tp+d$fp+d$fn))
d$same_species <- d$target_species == d$source_model
d$color <- ifelse(d$same_species, SAME, CROSS)
d <- d[order(match(d$target_species,SOURCES),match(d$source_model,SOURCES)),]
values <- as.matrix(d[c("pearson_r", "HiLearn_AUROC", "micro_jaccard")])
stopifnot(all(is.finite(values)), all(values>=0 & values<=1))
mm <- function(x) unit(x,"mm")
txt <- function(label,x,y,size=7,bold=FALSE,just="centre") {
  grid.text(label,x=mm(x),y=mm(y),just=just,gp=gpar(fontfamily="Arial", fontsize=size,
    fontface=if(bold) "bold" else "plain", col="black"))
}
line <- function(x,y,col="black",pt=.4,lty="solid") {
  grid.lines(x=mm(x),y=mm(y),gp=gpar(col=col,lwd=pt/.75,lty=lty,lineend="round",linejoin="round"))
}
point <- function(x,y,pch,col,size=1.8) {
  grid.points(x=mm(x),y=mm(y),pch=pch,size=mm(size),gp=gpar(col=col,fill="white",lwd=.85/.75))
}
title <- function(t) switch(t,
  arabidopsis=expression(plain("Input DNA: ")*italic("A. thaliana")),
  maize=expression(plain("Input DNA: ")*italic("Z. mays")),
  wheat=expression(plain("Input DNA: ")*italic("T. aestivum")))
coords <- list()
draw_radar <- function(record=FALSE, panels=1:3, legend=TRUE) {
  stopifnot(all(panels %in% 1:3), !anyDuplicated(panels))
  if (legend) {
  line(c(48,54),c(137,137),SAME,1.3);txt("Same species",56,137,just="left")
  line(c(103,109),c(137,137),CROSS,1.3);txt("Cross species",111,137,just="left")
  for (i in seq_along(SOURCES)) {
    x <- c(12,68,123)[i];model <- SOURCES[i]
    line(c(x,x+8),rep(129,2),"black",1.05,LTY[[model]])
    point(x+4,129,PCH[[model]],"black")
    txt(NAMES[[model]],x+11,129,just="left")
  }
  }
  for (i in panels) {
    t <- SOURCES[i];cx <- c(31.5,91.5,151.5)[i];cy <- 77;radius <- 25
    panel <- d[d$target_species==t,]
    txt(letters[i],cx-27.5,118,8,TRUE,just="left")
    txt(title(t),cx,118,7)
    txt(paste0("n = ",format(panel$n_observations_r[1],big.mark=",",scientific=FALSE)),cx,112.5,6.5)
    theta <- c(90,-30,210)*pi/180
    project <- function(v) cbind(x=cx+radius*v*cos(theta),y=cy+radius*v*sin(theta))
    for (v in c(.25,.5,.75,1)) {
      xy <- project(rep(v,3));line(c(xy[,1],xy[1,1]),c(xy[,2],xy[1,2]),if(v==1) "#AAAAB2" else "#DDDEE3",if(v==1) .5 else .35)
    }
    endpoints <- project(rep(1,3))
    for (j in 1:3) line(c(cx,endpoints[j,1]),c(cy,endpoints[j,2]),"#C6C7CE",.35)
    txt("Pooled Pearson r",cx,106,6.7)
    txt("Micro Jaccard (P85)",cx-26,58.5,6.4,just="left")
    txt("Micro AUROC",cx+26,58.5,6.4,just="right")
    # The grid interval is documented below; endpoint labels keep curves unobscured.
    txt("1",cx-1.1,cy+radius,5.8,just="right")
    txt("1",endpoints[2,1]+1.3,endpoints[2,2],5.8,just="left")
    txt("1",endpoints[3,1]-1.3,endpoints[3,2],5.8,just="right")
    txt("0",cx-1.2,cy-1,5.8,just="right")
    # Draw the own-species baseline last. No polygon fill and no area aggregation.
    for (model in c(setdiff(SOURCES,t),t)) {
      row <- panel[panel$source_model==model,]
      scores <- unlist(row[c("pearson_r","HiLearn_AUROC","micro_jaccard")],use.names=FALSE)
      xy <- project(scores);col <- row$color
      line(c(xy[,1],xy[1,1]),c(xy[,2],xy[1,2]),col,1.15,LTY[[model]])
      point(xy[,1],xy[,2],PCH[[model]],col)
      if(record) coords[[length(coords)+1]] <<- data.frame(target_species=t,source_model=model,
        metric=c("pearson_r","HiLearn_AUROC","micro_jaccard"),value=scores,x_mm=xy[,1],y_mm=xy[,2],
        center_x_mm=cx,center_y_mm=cy,radius_mm=radius,color=col)
    }
    # A small table supports precise reading of close model profiles.
    left <- cx-27;cols <- cx+c(4,14.2,24.5)
    txt("Source model",left,49,6.3,TRUE,just="left")
    for (j in 1:3) txt(c("r","AUROC","J85")[j],cols[j],49,6.3,TRUE)
    line(c(left,cx+28),rep(46.7,2),"#D1D2D8",.4)
    for (j in 1:3) {
      row <- panel[panel$source_model==SOURCES[j],];y <- 42.8-(j-1)*6
      point(left+1,y,PCH[[SOURCES[j]]],row$color,size=1.5)
      txt(NAMES[[SOURCES[j]]],left+3.2,y,6.3,just="left")
      scores <- unlist(row[c("pearson_r","HiLearn_AUROC","micro_jaccard")],use.names=FALSE)
      for(k in 1:3) txt(sprintf("%.3f",scores[k]),cols[k],y,6.3)
    }
  }
}


FONT <- "Arial"; BODY <- 7; PANEL <- 8
EXPERIMENT <- SAME; PREDICTION <- CROSS
FIG <- file.path(OUTPUT_DIR,"figures"); VAL <- file.path(OUTPUT_DIR,"validation")
track_map <- list()
rect <- function(x,y,w,h,color) grid.rect(x=mm(x),y=mm(y),width=mm(w),height=mm(h),just=c("left","bottom"),gp=gpar(fill=color,col=NA))
species_text <- title
regions <- read.csv(file.path(DATA_DIR,"track_selected.csv"),stringsAsFactors=FALSE)
limits <- read.csv(file.path(DATA_DIR,"track_y_limits.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(correlations)==9,nrow(peaks)==9,nrow(regions)==3,all(regions$display_rank==1))
inputs <- list.files(DATA_DIR,recursive=TRUE,full.names=TRUE)
before <- tools::md5sum(inputs)
tracks <- setNames(lapply(regions$candidate_id,function(id) {
 d <- read.csv(gzfile(file.path(DATA_DIR,"regions",paste0(id,".tracks.csv.gz"))),stringsAsFactors=FALSE)
 d$valid <- as.character(d$valid) %in% c("1","TRUE","True","true")
 stopifnot(all(d$input_end0-d$input_start0==1024),all(d$center_start0==d$input_start0+448),all(d$center_end0==d$input_start0+576))
 d
}),regions$candidate_id)
draw_tracks <- function(record=FALSE, panels=1:3) {
 stopifnot(all(panels %in% 1:3), !anyDuplicated(panels))
 for(i in panels) {
  target <- SOURCES[i];r <- regions[regions$target_species==target,,drop=FALSE];stopifnot(nrow(r)==1)
  d <- tracks[[r$candidate_id]];strand <- r$primary_strand;header <- c(86,58,30)[i]
  txt(letters[i+3],2,header,PANEL,TRUE,just="left")
  txt(species_text(target),8,header,just="left")
  txt(sprintf("Own-model r = %.3f",r$r_self),80,header,just="left")
  coord <- sprintf("Chr%s:%s-%s | %s",r$chrom,r$region_start0,r$region_end0,if(strand=="plus") "+" else "-")
  txt(coord,181,header,just="right",size=6.5)
  models <- c(target,r$cross_model1,r$cross_model2)
  for(j in 1:3) {
   experimental <- j==1;model <- models[j];track_color <- if(experimental) EXPERIMENT else PREDICTION;bottom <- header-c(8,14,20)[j];height <- 5;left <- 43;width <- 138
   column <- if(experimental) paste0("mean_cpm_128_",strand) else paste0("prediction_",model,"_",strand)
   values <- d[[column]];if(experimental) values <- log2(1+values)
   ymax <- if(experimental) limits$experimental_ymax[limits$target_species==target] else 1
   ok <- d$valid & is.finite(values);stopifnot(all(values[ok]<=ymax+1e-10))
   x <- left+(d$center_start0[ok]-r$region_start0)/(r$region_end0-r$region_start0)*width
   w <- (d$center_end0[ok]-d$center_start0[ok])/(r$region_end0-r$region_start0)*width
   rect(x,rep(bottom,length(x)),w,values[ok]/ymax*height,track_color)
   line(c(left,left,left+width),c(bottom+height,bottom,bottom),"#666666",.3)
   txt(format(ymax,trim=TRUE),left-1.2,bottom+height,size=6.5,just=c("right","top"))
   txt("0",left-1.2,bottom,size=6.5,just=c("right","bottom"))
   value_r <- NA_real_
   if(experimental) label <- "GRO-seq [log2(1+mean CPM)]" else {
    value_r <- cor(log2(1+d[[paste0("mean_cpm_1024_",strand)]][ok]),log2(1+d[[column]][ok]))
    expected <- if(j==2) r$r_model1 else r$r_model2;stopifnot(abs(value_r-expected)<1e-12)
    label <- sprintf("%s model, r = %.3f",NAMES[[model]],value_r)
   }
   txt(label,2,bottom+height/2,just="left")
   if(record) track_map[[length(track_map)+1L]] <<- data.frame(panel=letters[i+3],candidate_id=r$candidate_id,target=target,strand=strand,track=if(experimental) "experimental" else model,source_column=column,n_valid=sum(ok),n_invalid=sum(!d$valid),r=value_r,ymax=ymax,height_mm=height,width_mm=width,color=track_color)
  }
  ticks <- seq(r$region_start0,r$region_end0,length.out=3);xp <- c(43,112,181);base <- header-20
  for(k in 1:3) {line(rep(xp[k],2),c(base,base-.55),"#666666",.3);txt(sprintf("%.3f",ticks[k]/1e6),xp[k],header-22.1,size=6.5,just=c("left","centre","right")[k])}
  txt("Mb",35.5,header-22.1,size=6.5,just="right")
 }
}
set_exact_pdf_mediabox <- function(path,width_mm,height_mm) {
  original <- readBin(path,what="raw",n=file.info(path)$size)
  ascii <- original; values <- as.integer(ascii)
  ascii[(values<32L & !(values %in% c(9L,10L,13L))) | values>126L] <- as.raw(32L)
  txt <- rawToChar(ascii)
  stopifnot(nchar(txt,type="bytes")==length(original))
  hits <- gregexpr("startxref[[:space:]]+[0-9]+[[:space:]]+%%EOF",txt,perl=TRUE)[[1]]
  if(hits[[1]]<0) stop("Missing PDF startxref")
  k <- length(hits); last_hit <- hits[k]
  xref_text <- substr(txt,last_hit,last_hit+attr(hits,"match.length")[k]-1L)
  previous_xref <- as.numeric(regmatches(xref_text,regexpr("[0-9]+",xref_text,perl=TRUE)))
  if(substr(txt,previous_xref+1,previous_xref+4)!="xref") stop("Expected classic Cairo xref table")
  trailer_start <- regexpr("trailer[[:space:]]*<<",substr(txt,previous_xref+1,last_hit-1L),perl=TRUE)[[1]]
  if(trailer_start<0) stop("Missing PDF trailer")
  trailer <- substr(txt,previous_xref+trailer_start,last_hit-1L)
  if(grepl("/Encrypt",trailer,fixed=TRUE)) stop("Encrypted PDFs not supported")
  extract <- function(pattern,value=trailer) {
    found <- regmatches(value,regexpr(pattern,value,perl=TRUE))
    if(!length(found) || !nzchar(found)) stop("Missing PDF field: ",pattern)
    found
  }
  root <- extract("/Root[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]+R")
  size <- as.integer(sub("/Size[[:space:]]+","",extract("/Size[[:space:]]+[0-9]+")))
  info <- if(grepl("/Info",trailer,fixed=TRUE)) extract("/Info[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]+R") else ""
  identifier <- if(grepl("/ID",trailer,fixed=TRUE)) extract("/ID[[:space:]]*\\[[^]]*\\]") else ""
  starts <- gregexpr("(?m)^[0-9]+[[:blank:]]+[0-9]+[[:blank:]]+obj\\b",txt,perl=TRUE)[[1]]
  page_objects <- list()
  for(start in starts[starts>0]) {
    tail <- substr(txt,start,nchar(txt,type="bytes"))
    stop_at <- regexpr("endobj",tail,fixed=TRUE)[[1]]
    if(stop_at<0) next
    object <- substr(tail,1,stop_at+5L)
    if(grepl("/Type[[:space:]]*/Page\\b",object,perl=TRUE)) page_objects[[length(page_objects)+1L]] <- object
  }
  if(!length(page_objects)) stop("No PDF Page objects found")
  updates <- ""; entries <- list()
  box_pattern <- "/MediaBox[[:space:]]*\\[[^]]*\\]"
  exact_box <- sprintf("/MediaBox [ 0 0 %.12f %.12f ]",width_mm*72/25.4,height_mm*72/25.4)
  for(page in page_objects) {
    header <- strsplit(extract("^[0-9]+[[:blank:]]+[0-9]+[[:blank:]]+obj",page),"[[:blank:]]+")[[1]]
    number <- as.integer(header[1]); generation <- as.integer(header[2])
    stopifnot(size>number)
    extract(box_pattern,page)
    offset <- length(original)+nchar(updates,type="bytes")+1
    updates <- paste0(updates,"\n",sub(box_pattern,exact_box,page,perl=TRUE),"\n")
    entries[[length(entries)+1L]] <- data.frame(number=number,generation=generation,offset=offset)
  }
  entries <- do.call(rbind,entries); entries <- entries[order(entries$number),,drop=FALSE]
  stopifnot(!anyDuplicated(entries$number))
  new_xref <- length(original)+nchar(updates,type="bytes")
  xref <- "xref\n0 1\n0000000000 65535 f \n"
  for(i in seq_len(nrow(entries))) xref <- paste0(xref,entries$number[i]," 1\n",sprintf("%010d %05d n \n",entries$offset[i],entries$generation[i]))
  increment <- paste0(updates,xref,"trailer\n<< /Size ",size," ",root," ",info," ",identifier,
    " /Prev ",sprintf("%.0f",previous_xref)," >>\nstartxref\n",sprintf("%.0f",new_xref),"\n%%EOF\n")
  corrected <- c(original,charToRaw(increment))
  stopifnot(identical(corrected[seq_along(original)],original))
  writeBin(corrected,path)
  invisible(length(page_objects))
}


wrap_caption <- function(text,width_mm=183,size=8) {
 words <- strsplit(text,"[[:space:]]+")[[1]];lines <- character();current <- ""
 for(word in words) {
  candidate <- if(nzchar(current)) paste(current,word) else word
  w <- grid::convertWidth(grid::grobWidth(grid::textGrob(candidate,gp=grid::gpar(fontfamily=FONT,fontsize=size))),"mm",valueOnly=TRUE)
  if(w>width_mm && nzchar(current)) {lines <- c(lines,current);current <- word} else current <- candidate
 }
 c(lines,current)
}

caption_text <- readLines(file.path(ROOT,"figure_legend_submission.txt"),encoding="UTF-8",warn=FALSE)
caption_title <- caption_text[1]
caption_body <- paste(caption_text[nzchar(caption_text)][-1],collapse=" ")
caption_layout <- NULL
draw_artwork <- function(record=FALSE) {
  pushViewport(viewport(x=mm(0),y=mm(66),width=mm(183),height=mm(144),just=c("left","bottom"),clip="off"))
  draw_radar(record)
  popViewport()
  draw_tracks(record)
}
draw_artwork_page <- function(record=FALSE) {
  grid.newpage()
  pushViewport(viewport(x=mm(0),y=mm(0),width=mm(183),height=mm(210),just=c("left","bottom"),clip="off"))
  draw_artwork(record)
  popViewport()
}
draw_caption_page <- function(record=FALSE) {
  grid.newpage()
  pushViewport(viewport(x=mm(13.5),y=mm(77),width=mm(183),height=mm(210),just=c("left","bottom"),clip="off"))
  draw_artwork(record)
  popViewport()
  txt(caption_title,13.5,70.5,size=9,bold=TRUE,just=c("left","top"))
  lines <- wrap_caption(caption_body)
  ys <- 65.7-(seq_along(lines)-1)*10*25.4/72
  stopifnot(tail(ys,1)-8*25.4/72>=8)
  for(i in seq_along(lines)) txt(lines[i],13.5,ys[i],size=8,just=c("left","top"))
  if(record) caption_layout <<- data.frame(line=seq_along(lines),text=lines,y_mm=ys,font_pt=8,leading_pt=10)
}
render_composite <- function() {
if (!isTRUE(capabilities("cairo"))) stop("This figure requires R with Cairo graphics support")
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dir.create(VAL,recursive=TRUE,showWarnings=FALSE)
coords <<- list(); track_map <<- list(); caption_layout <<- NULL
write.csv(d,file.path(VAL,"radar_values.csv"),row.names=FALSE)
stem <- file.path(FIG,"GROseq_radar_tracks_with_legend_A4")
cairo_pdf(paste0(stem,".pdf"),width=210/25.4,height=297/25.4,family=FONT,bg="white")
draw_caption_page(TRUE);dev.off()
set_exact_pdf_mediabox(paste0(stem,".pdf"),210,297)
png(file.path(FIG,"GROseq_radar_tracks_with_legend_preview.png"),width=210,height=297,units="mm",res=200,type="cairo",bg="white")
draw_caption_page();dev.off()
art <- file.path(FIG,"GROseq_radar_tracks_artwork_183x210mm")
cairo_pdf(paste0(art,".pdf"),width=183/25.4,height=210/25.4,family=FONT,bg="white")
draw_artwork_page();dev.off()
set_exact_pdf_mediabox(paste0(art,".pdf"),183,210)
png(paste0(art,".png"),width=183,height=210,units="mm",res=600,type="cairo",bg="white")
draw_artwork_page();dev.off()
vertices <- do.call(rbind,coords)
vertices$artwork_x_mm <- vertices$x_mm
vertices$artwork_y_mm <- vertices$y_mm+66
write.csv(vertices,file.path(VAL,"radar_vertices.csv"),row.names=FALSE)
write.csv(do.call(rbind,track_map),file.path(VAL,"tracks_rendered.csv"),row.names=FALSE,na="")
write.csv(caption_layout,file.path(VAL,"caption_layout.csv"),row.names=FALSE)
stopifnot(identical(before,tools::md5sum(inputs)))
write.csv(data.frame(path=substring(names(before),nchar(DATA_DIR)+2L),md5=as.character(before),unchanged=TRUE),file.path(VAL,"input_integrity.csv"),row.names=FALSE)
message("Created vector A4 figure with caption and separate 183 x 210 mm artwork.")

}
if (!DEFINE_ONLY) render_composite()
