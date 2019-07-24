#' @useDynLib rLVIS
#' @import Rcpp methods
#' @export
processFloWave <- function(input, output) {
  processFloWave2(input, output)
}
