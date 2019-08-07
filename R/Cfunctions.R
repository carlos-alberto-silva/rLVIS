#' Extract metrics from LVIS Level 1 waves data
#'
#' This function will perform multiple tasks to allow for extracting
#' metrics from the full waveform data. It will first denoise and
#' fit gaussians to the waveform, find the ground values
#'
#' @useDynLib rLVIS
#' @import Rcpp methods
#' @export
level1ExtractMetrics <- function(input, output) {
  processFloWave2(input, output)
}
