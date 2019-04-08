#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'@usage readLevel1b(level1bpath)
#'
#'@param level1bpath file path pointing to LVIS level1b data (H5 format)
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'@examples
#'
#'# LVIS level1b file path
#'level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
#'unzip(level1b_filepath_zip, exdir = tempdir())
#'level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")
#'
#'# Reading LVIS level1b file
#'level1b<-readLevel1b(level1bpath=level1b_filepath)
#'
#'
#'@importFrom h5 h5file
#'@export
readLevel1b<-function(level1bpath) {
  Level1b<- h5::h5file(level1bpath, 'a')
  return(Level1b)
}



#'Read LVIS Level2 data
#'
#'@description This function reads LVIS level2 data
#'
#'@usage readLevel2(level2path, spdf=TRUE, glatlon=TRUE)
#'
#'@param level2path file path pointing to LVIS level2 data (txt format)
#'@param spdf if true, outups an object of class \code{SpatialPointsDataFrame}
#'@param glatlon if true, GLON and GLAT will be used for creating the \code{SpatialPointsDataFrame}.
#'If false, TLON and TLAT will be used instead
#'@return An object of class \code{SpatialPointsDataFrame} or \code{data.table};
#'@author Carlos Alberto Silva.
#'@examples
#'
#'# LVIS level2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'# Reading LVIS level1b file
#'level2_spdf<-readLevel2(level2path=level2_filepath, spdf=TRUE, glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'            colorPalette = c("blue","green","yellow","red"))
#'
#'
#'@export
readLevel2<-function(level2path, spdf=TRUE, glatlon=TRUE) {

  l2file<-utils::read.table(level2path,sep="")
  colnames(l2file)<-c("LFID","SHOTNUMBER","TIME","GLON","GLAT","ZG","TLON","TLAT","ZT","RH10","RH15","RH20","RH25",
                      "RH30","RH35","RH40","RH45","RH50","RH55","RH60","RH65","RH70","RH75","RH80","RH85",
                      "RH90","RH95","RH96","RH97","RH98","RH99","RH100","AZIMUTH","INCIDENTANGLE","RANGE","FLAG1",
                      "FLAG2","FLAG3")
  if (spdf==TRUE){

    if (glatlon==TRUE) {
      Level2<-sp::SpatialPointsDataFrame(l2file[,c("GLON", "GLAT")],data=l2file)
      sp::proj4string(Level2) <- sp::CRS("+proj=longlat +datum=WGS84")

    } else {

      Level2<-sp::SpatialPointsDataFrame(l2file[,c("TLON", "TLAT")],data=l2file)
      sp::proj4string(Level2) <- CRS("+proj=longlat +datum=WGS84")

    }
  } else {

    Level2<-data.table::data.table(l2file)

  }
  return(Level2)
}



#'Display LVIS Waveform
#'
#'@description This function plots LVIS level1b and level 2 data
#'
#'@usage plotWaveform(level1b,level2,shotnum,plot2=TRUE,...)
#'
#'@param level1b h5file; S4 object of class H5File
#'@param level2 dataframe containing LVIS level 2 data
#'@param shotnum LVIS shot number to display
#'@param plot2 if TRUE, plot both Level1b and Level2. If FALSE, plot only Level1b data
#'@param ... passing arguments on to the plot function
#'@return Returns 2-D scatterplot of the LVIS waveform
#'@author Carlos Alberto Silva.
#'@examples
#'
#'# LVIS level1b file path
#'level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
#'unzip(level1b_filepath_zip, exdir = tempdir())
#'level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")
#'
#'# Reading LVIS level1b file
#'level1b<-readLevel1b(level1bpath=level1b_filepath)
#'
#'# LVIS level2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'# Reading LVIS level1b file
#'level2_spdf<-readLevel2(level2path=level2_filepath, spdf=TRUE, glatlon=TRUE)
#'
#'#'Plotting LVIS waveforms
#'plotWaveform(level1b=level1b,level2=level2_spdf,
#'              shotnum=10964985,plot2=TRUE,xlab="Relative amplitude (%)", ylab="Height (m)")
#'
#'@importFrom grDevices colorRamp colors rgb
#'@importFrom graphics abline grid legend par plot polygon
#'@imortFrom utils read.table
#'@export
plotWaveform<-function(level1b,level2,shotnum=10964985,plot2=TRUE,...) {

  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }

  par(cex.axis=2)
  all_shotnums<-level1b['SHOTNUMBER'][]
  wave_idx<-which(all_shotnums==shotnum)
  wave_idx[wave_idx==TRUE]
  waveform<-level1b['RXWAVE'][wave_idx][1:1024]

  waveform<-(waveform/max(waveform))*100

  Z0<-level1b['Z0'][wave_idx]
  Z1023<-level1b['Z1023'][wave_idx]
  zstretch=rev(seq(Z1023,Z0,(Z0-Z1023)/1023))
  mytxtpathsub<-subset(level2@data,level2@data$SHOTNUMBER==shotnum)

  ZG<-mytxtpathsub[,6]
  ZT<-mytxtpathsub[,9]

  zstretch<-zstretch - ZG
  ZT<-ZT-ZG
  ZG<-0

  RH10<-mytxtpathsub[,10] #+ ZG
  RH25<-mytxtpathsub[,13] #+ ZG
  RH50<-mytxtpathsub[,18] #+ ZG
  RH75<-mytxtpathsub[,23] #+ ZG
  RH98<-mytxtpathsub[,30] #+ ZG
  RH100<-mytxtpathsub[,32] #+ ZG
  zmin<-ZG-(ZT-RH10)/6
  zmax=ZT + (ZT - RH10) /18
  x=zstretch>=zmin
  y=zstretch<=zmax
  z<-x==y
  waveform_crop=NULL
  zstrech_crop=NULL

  for ( i in 1:length(z)){
    if ( z[i]==TRUE) {
      waveform_crop[i]<-level1b['RXWAVE'][wave_idx][i]
      zstrech_crop[i]<-zstretch[i]
    }
  }

  #windows()
  if(plot2==TRUE) {
    par(mfrow=c(1,2))
    par(cex.axis=1.5)
    plot(waveform,zstretch, type="l", lwd=2, col="forestgreen")
    grid()
    polygon(c(waveform,min(zstretch),min(zstretch)),c(zstretch,min(zstretch),max(zstretch)),col="forestgreen")
    abline(h=ZG, lwd=2, col="blue")
    abline(h=RH10, lwd=2, col="green")
    abline(h=RH25, lwd=2, col="yellow")
    abline(h=RH50, lwd=2, col="gray")
    abline(h=RH75, lwd=2, col="orange")
    abline(h=RH98, lwd=2, lty=2,col="black")
    abline(h=RH100, lwd=2, col="red")

  }

  #browser()
  waveform_crop<-(as.numeric(waveform_crop)/as.numeric(summary(waveform_crop))[6])*100
  waveform_crop<- waveform_crop-as.numeric(summary(waveform_crop))[1]

  #polygon(c(0,as.numeric(summary(waveform_crop))[1],waveform_crop,as.numeric(summary(waveform_crop))[1],0,0),
  #        c(as.numeric(summary(zstrech_crop))[1],as.numeric(summary(zstrech_crop))[1],zstrech_crop,as.numeric(summary(zstrech_crop))[6],as.numeric(summary(zstrech_crop))[6],
  #          as.numeric(summary(zstrech_crop))[1]),col="forestgreen")

  xstart<-waveform_crop[which(zstrech_crop==min(zstrech_crop, na.rm=T))]
  xend<-waveform_crop[which(zstrech_crop==max(zstrech_crop, na.rm=T))]

  xl<-c(0,0,xstart,rev(waveform_crop),xend,0)
  yl<-c(max(zstrech_crop, na.rm=T),min(zstrech_crop, na.rm=T),min(zstrech_crop, na.rm=T),rev(zstrech_crop),max(zstrech_crop, na.rm=T),max(zstrech_crop, na.rm=T))

  par(cex.axis=1.5)
  plot(xl,yl, lwd=2, col="forestgreen", type="l",...)#,xlim=c(0,100))# ylim=c(min(zstretch),max(zstretch)))
  grid()
  #points(xl,yl)
  polygon(xl,yl,col="forestgreen")

  abline(h=ZG, lwd=2, col="blue")
  abline(h=RH10, lwd=2, col="green")
  abline(h=RH25, lwd=2, col="yellow")
  abline(h=RH50, lwd=2, col="gray")
  abline(h=RH75, lwd=2, col="orange")
  abline(h=RH98, lwd=2, lty=2,col="black")
  abline(h=RH100, lwd=2, col="red")
  par(xpd=TRUE)
  add_legend("topright",horiz=T,legend=c("ZG","RH10","RH25","RH50","RH75","RH98","RH100"),
             col=c("blue","green","yellow","gray","orange","black","red"),lwd=2,lty=c(rep(1,5),2,1),bty="n",cex=0.8)
  par(xpd=FALSE)
  #return(data.frame(SHOTNUMBER=rep(shotnum,1024),waveform,zstretch,min(yl, na.rm=T),max(yl, na.rm=T)))
}


#'LVIS Level2 data visualization in 2-D
#'
#'@description This function plots LVIS Level2 in 2-D
#'
#'@usage plotLevel2(level2_spdf, color, colorPalette,...)
#'
#'@param level2_spdf LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param color  The field name used to color the points.  Default is RH100
#'@param colorPalette A vector defining the color  palette.
#'@param ... passing arguments on to the plot function
#'@return A 2-D figure of LVIS Level2 data;
#'@author Carlos Alberto Silva.
#'@examples
#'
#'# LVIS level2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'# Reading LVIS level1b file
#'level2_spdf<-readLevel2(level2path=level2_filepath, spdf=TRUE, glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'head(level2_spdf@data)
#'#plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'#           colorPalette = c("blue","green","yellow","red"),
#'#           axes=T)
#'grid()
#'@export
plotLevel2 = function(level2_spdf, color = "RH100", colorPalette = c("blue","green","yellow","red"),...)
{
  v <- (level2_spdf@data[,color] - min(level2_spdf@data[,color]))/diff(range(level2_spdf@data[,color]))
  x <- colorRamp(colorPalette)(v)
  col<-rgb(x[,1], x[,2], x[,3], maxColorValue = 255)

  sp::plot(level2_spdf, col=col,...)
}


#'Statistics of LVIS Level2-derived Metrics
#'
#'@description Computes a Series of Statistics from LVIS Level2-derived Metrics
#'
#'@usage L2Stats(level2_spdf, func, id)
#'
#'@param level2_spdf LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param func the function to be applied to each cell
#'@param id a vector contatining the id for each LVIS Level2 observation. Defaut is NULL
#'@return Returns a \code{data.table} object containting Statistics of LVIS Level2-derived Metrics
#'@author Carlos Alberto Silva (This function has been adapted from the grid_metrics function in lidR package. All credits to Roussel et al. 2019).
#'@examples
#'
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'#' Polygons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_polygons.shp", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'           colorPalette = c("blue","green","yellow","red"))
#'
#'#' Reading Polygons
#'library(rgdal)
#'plots<-readOGR(polygons_filepath)
#'proj4string(plots) <- CRS("+proj=longlat +datum=WGS84")
#'plot(plots, add=TRUE, border="black", lwd=2)
#'
#'#'Clipping LVIS Level2 data
#'level2_spdf_sub<-clipLevel2(level2_spdf=level2_spdf,polygon_spdf=plots)
#'
#'#' Define your own function
#'mySetOfMetrics = function(x)
#'{
#'metrics = list(
#'    min =min(x), # Min of x
#'    max = max(x), # Max of x
#'    mean = mean(x), # Mean of x
#'    sd = sd(x)# Sd of x
#'  )
#'  return(metrics)
#'}
#'
#'#'Computing single LVIS metrics
#'RH100max<-L2Stats(level2_spdf=level2_spdf,func=~max(RH100), id=NULL)
#'
#'#'Computing LVIS metrics by id
#'RH100metrics<-L2Stats(level2_spdf=level2_spdf_sub,func=~mySetOfMetrics(RH100),
#'                      id=level2_spdf_sub@data$CLIPID)
#'
#'@export
L2Stats = function(level2_spdf, func=~mySetOfMetrics(RH100), id = NULL)
{
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  level2_dt<-data.table::data.table(level2_spdf@data)
  func<- lazyeval::f_interp(func)
  call<- lazyeval::as_call(func)

  if ( is.null(id)) {
    metrics   <- level2_dt[, c(eval(call))]
    metrics<-data.table::data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics)<-paste0(call)[1]
    }

  } else {
    metrics   <- level2_dt[, c(eval(call)), by = id]
    if (ncol(metrics) < 3) {
      colnames(metrics)[2]<-paste0(call)[1]
    }

  }

  return(metrics)
}

#'Compute a Series of Grid Metrics
#'
#'@description This function computes a series of user-defined descriptive statistics for a LVIS dataset within
#'each grid cell
#'
#'@usage GridMetrics(level2_spdf, func, res)
#'
#'@param level2_spdf LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param func the function to be applied to each cell
#'@param res spatial resolution for the output raster layer
#'@return Returns raster layer (s) of selected LVIS metric (s)
#'@author Carlos Alberto Silva (This function has been adapted from the grid_metrics function in lidR package. All credits to Roussel et al. 2019).
#'@examples
#'
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'#' Reading LVIS level 2 file
#'level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Define your own function
#'mySetOfMetrics = function(x)
#'{
#'metrics = list(
#'    min =min(x), # Min of z
#'    max = max(x), # Max of z
#'    mean = mean(x), # Mean of z
#'    sd = sd(x)# Sd of z
#'  )
#'  return(metrics)
#'}
#'
#'#'Computing LVIS metrics
#'mlvis<-GridMetrics(level2_spdf=level2_spdf,func=~mySetOfMetrics(ZT), res=0.0005)
#'
#'Computing single LVIS metrics
#'maxRH100<-GridMetrics(level2_spdf=level2_spdf,func=~max(RH100), res=0.0005)
#'
#'Computing single LVIS metrics
#'ZGmean<-GridMetrics(level2_spdf=level2_spdf,func=~mean(ZG), res=0.0005)
#'rasterVis::plot3D(ZGmean, col="gray")
#'
#'@importFrom utils read.table
#'@export
GridMetrics = function(level2_spdf, func=~max(RH100), res = 0.0005)
{
  # this code has been adapted from the grid_metrics function in lidR package (Roussel et al. 2019)
  # https://github.com/Jean-Romain/lidR/blob/master/R/grid_metrics.r

  level2_dt<-data.table::data.table(level2_spdf@data)
  layout<-raster::raster(raster::extent(level2_spdf), res=res)
  func      <- lazyeval::f_interp(func)
  call      <- lazyeval::as_call(func)
  cells     <- raster::cellFromXY(layout, sp::coordinates(level2_spdf))
  metrics   <- level2_dt[, c(eval(call)), by = cells]
  xy_coords <- raster::xyFromCell(layout, metrics[[1]])
  metrics[, cells := NULL]
  output <- sp::SpatialPixelsDataFrame(xy_coords, metrics, proj4string = level2_spdf@proj4string)
  names(output) <- names(metrics)
  if (length(names(metrics)) > 1 ) {output<-raster::brick(output)} else {output<-raster::raster(output)}
  return(output)
}

#'Clip LVIS Level2 data
#'
#'@description Clip LVIS Level2 data within a given geometry
#'
#'@usage clipLevel2(level2_spdf, polygon_spdf)
#'
#'@param level2_spdf h5file; S4 object of class H5File
#'@param polygon_spdf dataframe containing LVIS level 2 data
#'@return Returns An object of class \code{SpatialPoligonDataFrame} ; subset of LVIS Level2 data
#'@author Carlos Alberto Silva.
#'@examples
#'
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.txt")
#'
#'#' Polgons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_polygons.shp", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100", colorPalette = c("blue","green","yellow","red"))
#'
#'#' Reading Polygons
#'library(rgdal)
#'plots<-readOGR(polygons_filepath)
#'proj4string(plots) <- CRS("+proj=longlat +datum=WGS84")
#'plot(plots, add=T, border="black", lwd=2)
#'
#'#'Clipping LVIS Level2 data
#'level2_spdf_sub<-clipLevel2(level2_spdf=level2_spdf,polygon_spdf=plots)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf_sub, color = "RH100", colorPalette = c("blue","green","yellow","red"))
#'
#'@export
clipLevel2<-function(level2_spdf, polygon_spdf){
  CLIPID<-sp::over(level2_spdf,polygon_spdf)
  level2_spdf@data<-cbind(level2_spdf@data,CLIPID=CLIPID[,1])
  level2_spdf<-level2_spdf[rownames(CLIPID)[!is.na(CLIPID)],]
  return(level2_spdf)
}
