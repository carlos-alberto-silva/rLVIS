# Download GSL 2.4; hdf5; libgeotiff;

# Download gdal-2.2.0 from rwinlib
VERSION <- commandArgs(TRUE)
if(!file.exists(sprintf("../windows/gdal2-%s/include/gdal/gdal.h", VERSION))){
  download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/gsl-2.4/include/gsl/gsl_blas.h")){
  cat("Installing GSL...\n")
  download.file("https://github.com/rwinlib/gsl/archive/v2.4.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/mingw64-libhdf5-dev-1.8.20/include/hdf5.h")){
  cat("Installing HDF5...\n")
  download.file("https://github.com/caiohamamura/mingw64-libhdf5-dev/archive/v1.8.20.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

if(!file.exists("../windows/libgeotiff-1.4.3/geo_config.h")){
  cat("Installing libgeotiff...\n")
  download.file("https://github.com/OSGeo/libgeotiff/releases/download/1.4.3/libgeotiff-1.4.3.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
  file.rename("../windows/libgeotiff-1.4.3/geo_config.h.vc", "../windows/libgeotiff-1.4.3/geo_config.h")
}


if(!file.exists("../windows/libtiff-4.0.9/include/tiffio.h")){
  cat("Installing libtiff...\n")
  download.file("https://github.com/rwinlib/libtiff/archive/v4.0.9.zip", "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

source("../tools/deps.R")
