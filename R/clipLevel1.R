#'Clip LVIS Level1 data
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'@param level1_waveform h5file; S4 object of class H5File
#'@param output path where to save the new h5file
#'@param xleft numeric. left x coordinates of rectangles.
#'@param xright numeric. right x coordinates of rectangles.
#'@param ybottom numeric. bottom y coordinates of rectangles.
#'@param ytop numeric. top y coordinates of rectangles.
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@examples
#'
#'#' LVIS level 2 file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'# Rectangle
#'xleft = 9.35986
#'xright = 9.35988
#'ybottom = 0.5786
#'ytop = 0.5790
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1_filepath)
#'
#'output = tempfile(fileext="h5")
#'
#'clipped_waveform = clipLevel1(level1_waveform, output, xleft, xright, ybottom, ytop)
#'
#'@export
clipLevel1 = function(level1_waveform, output, xleft, xright, ybottom, ytop){

  spData = getSpatialData(level1_waveform)

  # xleft ybottom xright ytop
  mask =
    spData$lon0 >= xleft &
    spData$lon0 <= xright &
    spData$lat0 >= ybottom &
    spData$lat0 <= ytop &
    spData$lon1023 >= xleft &
    spData$lon1023 <= xright &
    spData$lat1023 >= ybottom &
    spData$lat1023 <= ytop

  mask = (1:length(spData$lon0))[mask]
  newFile = clipByMask(level1_waveform,
                       output,
                       mask)

  return (newFile)
}

#'Clip LVIS Level1 data by geometry
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'@param level1_waveform h5file; S4 object of class H5File
#'@param output path where to save the new h5file
#'@param polygon_spdf SpatialDataFrame. A polygon dataset for clipping the waveform
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@examples
#'
#'#' LVIS level 2 file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1_filepath)
#'
#'# Polgons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_clip_polygon.shp", package="rLVIS")
#'
#'# Reading LVIS level 2 file
#'polygon_spdf<-raster::shapefile(polygons_filepath)
#'
#'output = tempfile(fileext="h5")
#'
#'clipped_waveform = clipLevel1Geometry(level1_waveform, output, polygon_spdf)
#'
#'@export
clipLevel1Geometry = function(level1_waveform, output, polygon_spdf) {
  spData = getSpatialData(level1_waveform)

  points = sp::SpatialPointsDataFrame(coords=matrix(c(spData$lon0, spData$lat0), ncol=2),
                                      data=data.frame(id=1:length(spData$lon0)), proj4string = polygon_spdf@proj4string)
  pts = raster::intersect(points, polygon_spdf)
  mask = as.integer(pts@data$id)

  newFile = clipByMask(level1_waveform,
                       output,
                       mask)

  return (newFile)
}


getSpatialData = function(level1_waveform) {
  lat0 = level1_waveform["LAT0"][]
  lat1023 = level1_waveform["LAT1023"][]
  lon0 = level1_waveform["LON0"][]
  lon1023 = level1_waveform["LON1023"][]



  return (list(lat0=lat0, lon0=lon0, lat1023=lat1023, lon1023=lon1023))
}


clipByMask = function(level1_waveform, output, mask) {
  newFile = h5::h5file(output, mode="a")

  for (attr in h5::list.attributes(level1_waveform)) {
    h5::h5attr(newFile, attr) = h5::h5attr(level1_waveform, attr)
  }

  for (group in h5::list.groups(level1_waveform)) {
    g = newFile[group]
    for (dt in h5::list.datasets(level1_waveform[group], recursive = FALSE)) {
      newFile[dt] = level1_waveform[dt][]
    }
  }


  for (dataset in h5::list.datasets(level1_waveform, recursive = F)) {
    newFile[dataset] = level1_waveform[dataset][mask]
  }

  # Adjust ancillary extents if present
  if (h5::existsGroup(newFile, "ancillary_data")) {
    if (!h5::existsDataSet(newFile, "ancillary_data/Maximum Latitude")) {
      h5::createDataSet(newFile, "ancillary_data/Maximum Latitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Maximum Longitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Minimum Latitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Minimum Longitude", type="character", dimensions=1)
    }
    maxY = max(newFile["LAT0"][], newFile["LAT1023"][])
    maxX = max(newFile["LON0"][], newFile["LON1023"][])
    minY = min(newFile["LAT0"][], newFile["LAT1023"][])
    minX = min(newFile["LON0"][], newFile["LON1023"][])

    writeH5field(newFile["ancillary_data/Maximum Latitude"], sprintf("%.8f",maxY))
    writeH5field(newFile["ancillary_data/Maximum Longitude"], sprintf("%.8f",maxX))
    writeH5field(newFile["ancillary_data/Minimum Latitude"], sprintf("%.8f",minY))
    writeH5field(newFile["ancillary_data/Minimum Longitude"], sprintf("%.8f",minX))
  }

  return (newFile)
}


writeH5field = function(dataset, value) {
  ds = h5::selectDataSpace(dataset, 1)
  h5::writeDataSet(dataset, value, ds)
}
