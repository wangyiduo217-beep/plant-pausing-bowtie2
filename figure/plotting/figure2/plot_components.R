# Reuse the final figure's vector drawing functions; do not redraw from images.
THIS_DIR <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1],winslash="/")) else getwd()
})
GROSEQ_DEFINE_ONLY <- TRUE
source(file.path(THIS_DIR,"plot_radar_tracks_composite.R"),local=FALSE)
if (!isTRUE(capabilities("cairo"))) stop("This figure requires R with Cairo graphics support")
component_arg <- which(args == "--component")
if (length(component_arg)>1L || (length(component_arg) && component_arg==length(args))) stop("Use --component NAME once")
component <- if (length(component_arg)) args[component_arg+1L] else "all"
stopifnot(component %in% c("all","radar","tracks","legend",letters[1:6]))
COMPONENT_DIR <- file.path(OUTPUT_DIR,"components")
dir.create(COMPONENT_DIR,recursive=TRUE,showWarnings=FALSE)
rendered <- list()

render_crop <- function(name,width_mm,height_mm,x0_mm,y0_mm,draw) {
  draw_page <- function() {
    grid.newpage()
    pushViewport(viewport(x=mm(-x0_mm),y=mm(-y0_mm),width=mm(183),height=mm(144),
                          just=c("left","bottom"),clip="off"))
    draw()
    popViewport()
  }
  pdf <- file.path(COMPONENT_DIR,paste0(name,".pdf"))
  cairo_pdf(pdf,width=width_mm/25.4,height=height_mm/25.4,family=FONT,bg="white")
  draw_page();dev.off()
  set_exact_pdf_mediabox(pdf,width_mm,height_mm)
  png(file.path(COMPONENT_DIR,paste0(name,".png")),width=width_mm,height=height_mm,
      units="mm",res=200,type="cairo",bg="white")
  draw_page();dev.off()
  rendered[[length(rendered)+1L]] <<- data.frame(component=name,width_mm=width_mm,height_mm=height_mm,
      source_origin_x_mm=x0_mm,source_origin_y_mm=y0_mm,
      geometry="Original vector functions; translation only; no scaling or numerical changes")
}

if (component %in% c("all","radar")) {
  render_crop("radar_panels_a_c",183,144,0,0,function() draw_radar())
}
if (component %in% c("all","tracks")) {
  render_crop("track_panels_d_f",183,90,0,0,function() draw_tracks())
}
if (component %in% c("all","legend")) {
  render_crop("shared_radar_legend",183,22,0,122,function() draw_radar(panels=integer(0),legend=TRUE))
}
for (i in 1:3) {
  if (component %in% c("all",letters[i])) {
    local({
      panel_index <- i
      render_crop(paste0("panel_",letters[i],"_radar_",SOURCES[i]),63,98,(i-1)*60,26,
                  function() draw_radar(panels=panel_index,legend=FALSE))
    })
  }
  if (component %in% c("all",letters[i+3])) {
    local({
      panel_index <- i
      render_crop(paste0("panel_",letters[i+3],"_tracks_",SOURCES[i]),183,28,0,c(62,34,6)[i],
                  function() draw_tracks(panels=panel_index))
    })
  }
}
write.csv(do.call(rbind,rendered),file.path(COMPONENT_DIR,paste0("components_",component,".csv")),row.names=FALSE)
stopifnot(identical(before,tools::md5sum(inputs)))
message("Created vector figure components using the original drawing functions.")
