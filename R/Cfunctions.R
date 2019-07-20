#' @useDynLib rLVIS processFloWave2
#' @export
metrics <- function(input, output) {
  .Call("processFloWave2", input, output)
}
