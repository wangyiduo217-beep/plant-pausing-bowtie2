# Reproduce four R figure families. ROC/PR has a separate Python entry point.
run_publication_plots <- function() {
  own <- NULL
  for(fr in rev(sys.frames())) if(!is.null(fr$ofile) && basename(fr$ofile)=="run_all.R")
    {own <- dirname(normalizePath(fr$ofile,winslash="/",mustWork=TRUE));break}
  if(is.null(own)) {
    args <- commandArgs(FALSE);f <- sub("^--file=","",args[grepl("^--file=",args)])
    if(length(f))own <- dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE))
  }
  if(is.null(own))stop("Source the full path to figure/plotting/run_all.R.")
  entries <- c("correlations/run_plots.R","tracks/run_tracks.R",
               "regional_bars/run_regional_bars.R","deciles/run_decile_correlations.R")
  for(entry in entries) {
    message("Running: ",entry)
    saved_options <- options()
    tryCatch(source(file.path(own,entry),local=new.env(parent=globalenv()),encoding="UTF-8"),
             finally=options(saved_options))
  }
  message("Completed all four R figure families under figure/results/.")
}
run_publication_plots()
