# Draw the frozen strand-specific decile statistics; never reruns bootstrap.
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

module_dir <- find_module_dir("run_decile_correlations.R")
figure_dir <- normalizePath(file.path(module_dir,"../.."),winslash="/",mustWork=TRUE)
source(file.path(module_dir,"config.R"),local=environment(),encoding="UTF-8")
if(!DECILE_LAYOUT %in% c("grouped","sixpanel")) stop("DECILE_LAYOUT must be grouped or sixpanel")
if (!capabilities("cairo")) stop("Cairo graphics support is required for embedded text.")
statistics_path <- file.path(figure_dir,"data/deciles/group_statistics.csv")
output_dir <- file.path(figure_dir,"results/deciles")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source(file.path(module_dir,paste0("draw_",DECILE_LAYOUT,".R")),local=environment(),encoding="UTF-8")
