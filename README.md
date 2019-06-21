![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/fig_1.png)<br/>

![Github](https://img.shields.io/badge/CRAN-0.0.2-green.svg)
![Github](https://img.shields.io/badge/Github-0.0.2-green.svg)
[![Rdoc](http://www.rdocumentation.org/badges/version/rLVIS)](http://www.rdocumentation.org/packages/rLVIS)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![R_Forge](https://img.shields.io/badge/R_Forge-0.0.2-green.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rLVIS)

rLVIS: An R Package for NASA's Land, Vegetation, and Ice Sensor (LVIS) data processing and visualization

Set of tools for reading, processing and visualizing 
            LVIS (Land, Vegetation, and Ice Sensor) Level1b and Level2 data
            for forest inventory applications.
            
The rLVIS package provides functions to i) import Level1b and Level2 data, ii) visualing LVIS waveform, iii)
detect ground elevation from the waveform; iv) compute canopy metricfrom the LVIS waveform; v)  clip Level2 data within field plots; vi) to compute a series of statistics from Level2 at plot or as raster layers for the entire landscapy, vii) compare predicted version oberved forest attributes.

## Installation
```r
#The development version:
library(devtools)
devtools::install_github("carlos-alberto-silva/rLVIS")

#The CRAN version:
install.packages("rLVIS")
```    

## Getting Started
```r   
# Import labraries
library(rLIVS)
library(raster)
library(rgdal)
library(rasterVis)
library(rgl)
```

# Import LVIS data
```r
# LVIS level1b file path
level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
unzip(level1b_filepath_zip, exdir = tempdir())
level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")

#Reading LVIS level1b
level1b<-readLevel1b(level1bpath=level1b_filepath)
level1b

#' LVIS level 2 file path
level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
unzip(level2_filepath_zip, exdir = tempdir())
level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")

#' Reading LVIS level 2 file
level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)

#' Reading Polygons
polygons_filepath <- system.file("extdata", "LVIS_Mondah_polygons.shp", package="rLVIS")
plots<-readOGR(polygons_filepath)
proj4string(plots) <- CRS("+proj=longlat +datum=WGS84")
```

#' Plot LVIS Level2 data
```r
par(mfrow=c(1,2))
plotLevel2(level2_spdf=level2_spdf, color = "RH100",
           colorPalette = c("blue","green","yellow","red"),axes=TRUE,
           xlab="Lat", ylab="Lon")
plot(plots[3,], add=TRUE, border="black", lwd=2)

plotLevel2(level2_spdf=level2_spdf, color = "RH100",
           colorPalette = c("blue","green","yellow","red"),
           axes=TRUE, xlim=extent(plots[3,])[1:2],ylim=extent(plots[3,])[3:4],xlab="Lat", ylab="Lon")
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_1.png)


#'Clipping LVIS Level2 data
```r
level2_spdf_sub<-clipLevel2(level2_spdf=level2_spdf,polygon_spdf=plots)

#' Plot clipped LVIS Level2 data
par(mfrow=c(1,2))
plotLevel2(level2_spdf=level2_spdf, color = "RH100",
           colorPalette = c("blue","green","yellow","red"),axes=TRUE,
           xlab="Lat", ylab="Lon")
plot(plots, border="black", lwd=2, axes=T)
plotLevel2(level2_spdf=level2_spdf_sub, color = "RH100",
           colorPalette = c("blue","green","yellow","red"), add=T)
grid()
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_2.png)


#'Plotting LVIS waveforms
```r
plotWaveform(level1b=level1b,level2=level2_spdf,
             shotnum=10964985,plot2=TRUE,xlab="Relative amplitude (%)", ylab="Height (m)")
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_3.png)


#' Canopy and ground metrics within field plots
```r
#' Define your own function
mySetOfMetrics = function(x)
{
  metrics = list(
    min =min(x), # Min of z
    max = max(x), # Max of z
    mean = mean(x), # Mean of z
    sd = sd(x)# Sd of z
  )
  return(metrics)
}

#'Computing single LVIS metrics
RH100mean<-L2Stats(level2_spdf=level2_spdf,func=~mean(RH100), id=NULL)
head(RH100mean)                      
    ##    mean
    ## 34.5598
```
```r
#'Computing LVIS metrics by id
RH100metrics<-L2Stats(level2_spdf=level2_spdf_sub,func=~mySetOfMetrics(RH100),
                      id=level2_spdf_sub@data$CLIPID)
head(RH100metrics)                      

    ##   id  min   max   mean        sd
    ##   1 2.06 65.40 32.48820  9.996999
    ##   3 2.47 57.26 37.95028 12.054305
    ##   2 6.92 59.78 37.23889  5.176369
```
#' Canopy and ground metrics within as raster layers

```r
#' Computing serie of LVIS metrics
mlvis<-GridMetrics(level2_spdf=level2_spdf,func=~mySetOfMetrics(RH100), res=0.0005)
plot(mlvis)
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_4.png)

```r
#' Computing single LVIS metrics
maxRH100<-GridMetrics(level2_spdf=level2_spdf,func=~max(RH100), res=0.0005)
plot(maxRH100, xlab="UTM Easting", ylab="UTM Nothing")
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_5.png)

```r
plot3D(maxRH100, col="forestgreen")
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_7.PNG)

```r
#' Computing single LVIS metrics
ZGmean<-GridMetrics(level2_spdf=level2_spdf,func=~mean(ZG), res=0.0005)
plot(ZGmean, xlab="UTM Easting", ylab="UTM Nothing")
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_6.png)

```r
plot3D(ZGmean, col="gray", add=T)
aspect3d(1,1,0.1)
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_8.PNG)


#' Scatterplot of a 1:1 comparison

```r
#Importing libraries
library(raster)
library(rasterVis)
library(viridis)
library(gridExtra)

# Importing dataset
sf_agb1ha_path <- system.file("extdata", "sf_agb_1ha.tif", package="rLVIS")
lf_agb1ha_path <- system.file("extdata", "lf_agb_1ha.tif", package="rLVIS")

sf_agb<-raster(sf_agb1ha_path)
lf_agb<-raster(lf_agb1ha_path)

# Ploting AGB maps
s <- stack(sf_agb,lf_agb)
agb.maps<-levelplot(s,
          layout=c(1, 2),
          margin=FALSE,
          colorkey=list(
            space='right',
            labels=list(at=seq(0, 500, 50), font=4),
            axis.line=list(col='black'),
            width=1),
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(draw=TRUE),
          col.regions=viridis,
          at=seq(0, 500, len=101),
          names.attr=c("SF_AGB","LF_AGB"))

# Ploting Stats
colours<-viridis(10)
base_size=15
legend.position= c(0.85, 0.3)
stat.size=5
stats.position=c(100,400,50)
base_size=base_size
xlim=c(0,500)
legend.size=c(8,15,10,2)
ylim=c(0,500)
ylab="LF-derived AGB (Mg/ha)"
xlab="SF-derived AGB (Mg/ha)"
fit.line.col=c("black","gray")
title="SF vs LF lidar"
x_axis=TRUE
y_axis=TRUE

x11()
agb.comp<-plotStats(y=getValues(r1),
                   x=getValues(r2),
                   colours=colours,
                   legend.position= legend.position,
                   stat.size=stat.size,
                   stats.position=stats.position,
                   base_size=base_size,
                   xlim=xlim,
                   legend.size=legend.size,
                   ylim=ylim,
                   ylab=ylab,
                   xlab=xlab,
                   fit.line.col=fit.line.col,
                   title=title,
                   x_axis=x_axis,
                   y_axis=y_axis)

# Combining plots
grid.arrange(agb.maps,agb.comp$plotg, nrow = 1)
```
![](https://github.com/carlos-alberto-silva/rLVIS/blob/master/readme/Fig_9.PNG)
