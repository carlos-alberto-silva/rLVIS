#' @useDynLib rLVIS processFloWave2
#' @export
processFloWave <- function(input, output) {
  .Call("processFloWave2", input, output)
}
