if(!file.exists("./cmpfit-1.2/mpfit.h")) {
  cat("Downloading cmpfit-1.2...\n")
  download.file("https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz", "cmpfit.tar.gz", quiet=TRUE)
  untar("cmpfit.tar.gz", exdir=".")
  unlink("cmpfit.tar.gz")
}
