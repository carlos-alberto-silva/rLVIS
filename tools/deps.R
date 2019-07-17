if(!file.exists("../cmpfit/mpfit.h")) {
  download.file("https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz", "cmpfit.tar.gz", quit=TRUE)
  dir.create("../cmpfit", showWarnings = FALSE)
  untar("cmpfit.tar.gz", exdir="cmpfit")
  unlink("cmpfit.tar.gz")
}
