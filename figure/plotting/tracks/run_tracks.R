# Reuse the audited correlation/track composite engine with the latest track tables.
run_tracks <- function() {
  own <- NULL
  for(fr in rev(sys.frames())) if(!is.null(fr$ofile) && basename(fr$ofile)=="run_tracks.R")
    {own <- dirname(normalizePath(fr$ofile,winslash="/",mustWork=TRUE));break}
  if(is.null(own)) {
    args <- commandArgs(FALSE);f <- sub("^--file=","",args[grepl("^--file=",args)])
    if(length(f))own <- dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE))
  }
  if(is.null(own))stop("Source the full path to run_tracks.R.")
  names <- c("GROSEQ_BUNDLE_DIR","GROSEQ_CONFIG","GROSEQ_MIN_PREDICTION_SCORE")
  old <- Sys.getenv(names,unset=NA_character_)
  on.exit(for(i in seq_along(names)) {
    if(is.na(old[i]))Sys.unsetenv(names[i]) else do.call(Sys.setenv,setNames(list(old[i]),names[i]))
  },add=TRUE)
  engine <- normalizePath(file.path(own,"..","correlations"),winslash="/",mustWork=TRUE)
  Sys.setenv(GROSEQ_BUNDLE_DIR=engine,GROSEQ_CONFIG=file.path(own,"config.R"),GROSEQ_MIN_PREDICTION_SCORE="0")
  source(file.path(engine,"run_plots.R"),local=new.env(parent=globalenv()),encoding="UTF-8")
}
run_tracks()
