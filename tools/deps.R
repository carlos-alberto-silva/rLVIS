downloadDep = function(name, file, url) {
  if(!file.exists(paste0("./",name,"/",file))) {
    cat(paste0("Downloading ",name,"...\n"))
    download.file(url, "lib.tar.gz", quiet=TRUE)
    untar("lib.tar.gz", exdir=".")
    unlink("lib.tar.gz")
  }
}


downloadDepBitBucket = function(name, file, origName) {
  if(!file.exists(paste0("./",name,"/",file))) {
    cat(paste0("Downloading ",name,"...\n"))
    url = paste0("https://bitbucket.org/caiohamamura/",name,"/get/v0.1.2.zip")
    download.file(url, "lib.zip", quiet=TRUE)
    unzip("lib.zip", exdir=".")
    unlink("lib.zip")
    file.rename(origName, name)
  }
}

downloadDep("cmpfit-1.2",
            "mpfit.h",
            "https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz")
downloadDepBitBucket("gedisimulator",
                     "gediRat.c",
                     "caiohamamura-gedisimulator-4136c9ed7fba")
downloadDepBitBucket("tools",
                     "tools.c",
                     "caiohamamura-tools-b7ba3550791e")
downloadDepBitBucket("libclidar",
                     "libLasProcess.h",
                     "caiohamamura-libclidar-bf4cf5fe087a")
