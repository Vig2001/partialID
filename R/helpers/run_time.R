# Wall-clock run time of an R script, for the console and run_info.txt.


#' Time elapsed since start, as "1 h 23 min 4 s" (hours omitted if zero).
#'
#' @param start POSIXct, e.g. Sys.time() at the start of main().
format_elapsed <- function(start) {
  seconds <- round(as.numeric(difftime(Sys.time(), start, units = "secs")))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  if (hours > 0) {
    sprintf("%d h %d min %d s", hours, minutes, secs)
  } else {
    sprintf("%d min %d s", minutes, secs)
  }
}
