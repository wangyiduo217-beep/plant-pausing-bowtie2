# Run four fixed display thresholds without editing config.R or the original data.
locate_bundle <- function() {
  override <- Sys.getenv("GROSEQ_BUNDLE_DIR", "")
  if(nzchar(override)) return(normalizePath(override,winslash="/",mustWork=TRUE))
  frames <- sys.frames()
  for(i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if(!is.null(f) && file.exists(f)) {
      d <- dirname(normalizePath(f,winslash="/",mustWork=TRUE))
      if(file.exists(file.path(d,"run_plots.R"))) return(d)
    }
  }
  args <- commandArgs(trailingOnly=FALSE)
  f <- sub("^--file=","",args[grepl("^--file=",args)])
  if(length(f)) return(dirname(normalizePath(f[[1]],winslash="/",mustWork=TRUE)))
  if(file.exists("run_plots.R")) return(normalizePath(getwd(),winslash="/"))
  stop("Cannot find run_plots.R. Open GROseq_plotting.Rproj before running this file.")
}
run_comparison <- function() {
  bundle <- locate_bundle()
  settings <- new.env(parent=baseenv())
  sys.source(file.path(bundle,"config.R"),envir=settings)
  thresholds <- settings$COMPARISON_THRESHOLDS
  stopifnot(is.numeric(thresholds),all(is.finite(thresholds)),all(thresholds>=0 & thresholds<=1))
  old_threshold <- Sys.getenv("GROSEQ_MIN_PREDICTION_SCORE",unset=NA_character_)
  old_bundle <- Sys.getenv("GROSEQ_BUNDLE_DIR",unset=NA_character_)
  on.exit({
    if(is.na(old_threshold)) Sys.unsetenv("GROSEQ_MIN_PREDICTION_SCORE") else Sys.setenv(GROSEQ_MIN_PREDICTION_SCORE=old_threshold)
    if(is.na(old_bundle)) Sys.unsetenv("GROSEQ_BUNDLE_DIR") else Sys.setenv(GROSEQ_BUNDLE_DIR=old_bundle)
  },add=TRUE)
  Sys.setenv(GROSEQ_BUNDLE_DIR=bundle)
  tables <- list()
  for(i in seq_along(thresholds)) {
    Sys.setenv(GROSEQ_MIN_PREDICTION_SCORE=format(thresholds[i],scientific=FALSE,trim=TRUE,digits=12))
    env <- new.env(parent=globalenv())
    source(file.path(bundle,"run_plots.R"),local=env,encoding="UTF-8")
    tables[[i]] <- env$metrics
  }
  combined <- do.call(rbind,tables)
  write.csv(combined,file.path(bundle,"..","..","generated","correlations","threshold_comparison_statistics.csv"),row.names=FALSE)
  message("Finished fixed-threshold comparison: ",paste(thresholds,collapse=", "))
}
run_comparison()
