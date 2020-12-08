#Code to Process Landsat HLS data!
---
  # Transect data (ground truth) over several counties in Indiana that tell us cover crop presence/absence in addition to cover crop type over several fields. 
  #Specifically we have data for 5 counties: 
  # Posey County (Fall 2015, Fall 2016, Spring 2017)
  #Benton County (Fall 2014, Fall 2015, Fall 2016, Spring 2017)
  #Warren County (Fall )
  #White County
  #Gibson County (Fall 2015, Fall 2016, Spring 2017)
  
#Comparisons: 
# Building presence/absence models for cover crop detection
#Cover crop detection: NDVI -> winter integral only
#Cover crop detection: NDVI (ratio fall months to winter integral)
#Cover crop detection: VisNIR
#Cover crop detection: SWIR (including NDRI and NDTI)
#Cover crop detection: Thermal Data

#Harmonized Landast-Sentinel 2 data?: 
# https://hls.gsfc.nasa.gov/products-description/

pacman::p_load(raster, gdalUtils, terra, dplyr, rgdal, sp, rhdf5, rlist, randomForest, caTools, caret, maptools)

#  We are using the L30 Landsat-like product (and comparing it with Landsat) - spatial resolution 30 m, temporal resolution 5-day. Our sentinel tiles are: 16T, 17T, 16S, 17S (see here: https://hls.gsfc.nasa.gov/wp-content/uploads/2016/03/MGRS_GZD-1.png)
#Mean NA basic function for raster operations
mean_na <- function(x) {
  mean(x,na.rm=T)
}

#Options for terra and raster 
terraOptions(progress=10, memfrac=0.6)
rasterOptions(tmpdir="C:\\",tmptime = 24,progress="text",timer=TRUE,overwrite = T,chunksize=2e+08,maxmemory=1e+8)

#Files come as HDF
# We are using the L30 Landsat-like product (and comparing it with Landsat) - spatial resolution 30 m, temporal resolution 5-day. Our sentinel tiles are: 16T, 17T, 16S, 17S (see here: https://hls.gsfc.nasa.gov/wp-content/uploads/2016/03/MGRS_GZD-1.png)
#Function to project all to Gtiff

projHDF2GTiff = function(loc, hdfs, gtiffs, lyr, fromSRS, toSRS){ 
  if("gdalUtils" %in% rownames(installed.packages()) == FALSE){
    install.packages("gdalUtils", repos="http://r-forge.r-project.org")
    require(gdalUtils)
  } # install and load the gdalUtils package. 
  setwd(loc) # set the working directory to where the data is
  suppressWarnings(dir.create(paste(loc,"Projected", lyr,sep="/"))) # create a directory to store projected files
  for (i in 1:length(hdfs)){ 
    gdal_translate(hdfs[i],gtiffs[i],sd_index=lyr) # extract the specified HDF layer and save it as a Geotiff
    gdalwarp(gtiffs[i],paste(loc,"Projected",lyr, gtiffs[i],sep="/"),s_srs=fromSRS,t_srs=toSRS,srcnodata=-1000,dstnodata=NA,overwrite = T) # project geotiffs
    unlink(gtiffs[i]) # delete unprojected geotiffs to save space
  }
}

#x is tile, y is 'year'
Process_L30 <- function(x,y){
  #pass x to the "my loc" variable
  myloc=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, y, "L30", sep="/")
  print(myloc)
  setwd(myloc)
  hdfs1 = list.files(getwd(), pattern="hdf$")
  gtiffs1 = gsub("hdf","tif",hdfs1) #
  frm.srs = "+proj=utm +zone=16 +ellps=WGS84 + datum=WGS84 + units=m + no_defs" # original HDF SRS
  to.srs = "+proj=longlat +datum=WGS84 +no_defs" # desired GeoTIFF SRS
  # lyr is the HDF layer you want to extract. In this example it is "1" to 
  #Project to .tiff - do it with all bands - files will be named by L30 Subdataset number 
  print("Band1 to Tiff")
  #Band 1
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 1, fromSRS = frm.srs, toSRS = to.srs)
  #Band 2
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 2, fromSRS = frm.srs, toSRS = to.srs)
  #Band 3
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 3, fromSRS = frm.srs, toSRS = to.srs)
  #Band 4
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 4, fromSRS = frm.srs, toSRS = to.srs)
  #Band 5
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 5, fromSRS = frm.srs, toSRS = to.srs)
  #Band 6
  print("Band6 to Tiff")
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 6, fromSRS = frm.srs, toSRS = to.srs)
  #Band 7
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 7, fromSRS = frm.srs, toSRS = to.srs)
  #Band 9 
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 8, fromSRS = frm.srs, toSRS = to.srs)
  #Band 10 
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 9, fromSRS = frm.srs, toSRS = to.srs)
  #Band 11
  print("Band QA to Tiff")
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 10, fromSRS = frm.srs, toSRS = to.srs)
  #QA Band
  projHDF2GTiff(loc = myloc, hdfs = hdfs1, gtiffs = gtiffs1, lyr = 11, fromSRS = frm.srs, toSRS = to.srs)
  rm(myloc,hdfs1,gtiffs1,frm.srs,to.srs,srcnodata,dstnodata) # remove variables to save memory
#Now create raster stacks
  return()
  print("done")
  }
write_bandstacks <- function(x, d, c){
  #Let's create a raster stack for each 2015-2016 WY timestep...with all relevant bands (3,4,5,6,7,10,11)
  #D is the year and C is the year-1
  #Now for band 3....
  #Getting 2015 band 3
  print("band3")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/3", sep="/")
  setwd(myloc1)
  band3_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/3", sep="/")
  setwd(myloc2)
  band3_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band3_2016))
  print(length(band3_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band3_2015[c(35:45)],band3_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana",x, "2015-2016wy/Band_3.tif",sep="/")
  print(out)
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band4")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/4", sep="/")
  setwd(myloc1)
  band4_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/4", sep="/")
  setwd(myloc2)
  band4_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band4_2016))
  print(length(band4_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b42016 <- c(band4_2015[c(35:45)],band4_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_4.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band5")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/5", sep="/")
  setwd(myloc1)
  band5_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/5", sep="/")
  setwd(myloc2)
  band5_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band5_2016))
  print(length(band5_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band5_2015[c(35:45)],band5_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_5.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band6")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/6", sep="/")
  setwd(myloc1)
  band6_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/6", sep="/")
  setwd(myloc2)
  band6_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band6_2016))
  print(length(band6_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band6_2015[c(35:45)],band6_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_6.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band7")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/7", sep="/")
  setwd(myloc1)
  band7_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/7", sep="/")
  setwd(myloc2)
  band7_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band7_2016))
  print(length(band7_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band7_2015[c(35:45)],band7_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_7.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band9")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/9", sep="/")
  setwd(myloc1)
  band9_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/9", sep="/")
  setwd(myloc2)
  band9_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band9_2016))
  print(length(band9_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band9_2015[c(35:45)],band9_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_9.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band10")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/10", sep="/")
  setwd(myloc1)
  band9_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 3
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/10", sep="/")
  setwd(myloc2)
  band9_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band3_2016))
  print(length(band3_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band3_2015[c(35:45)],band3_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  scaled_tstack <- calc(tststack, function(x){x*0.0001})
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_10.tif",sep="/")
  writeRaster(scaled_tstack, out, overwrite=TRUE)
  
  print("band11 (QA)")
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, c, "L30/Projected/11", sep="/")
  setwd(myloc1)
  band3_2015 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  #Getting 2016 band 11
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, d, "L30/Projected/11", sep="/")
  setwd(myloc2)
  band3_2016 <- list.files(getwd(), pattern="tif$", full.names=TRUE)
  print(length(band3_2016))
  print(length(band3_2015))
  #Moving all files in the folloiwng list to
  print(paste("stacking", c, d, "to get water year", sep=" "))
  b32016 <- c(band3_2015[c(35:45)],band3_2016[c(1:34)])
  #All of band 3 for 2016 is stacked! 
  tststack <- raster::stack(b32016)
  #Changes the names to day of year
  names(tststack) <- paste("doy",substr(names(tststack), 20,22), sep="_")
  #Takes about 143 seconds
  out=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy/Band_11.tif",sep="/")
  writeRaster(tststack, out, overwrite=TRUE)
  gc()
  return()
  }
masking <- function(x,y){
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy",  sep="/")
  setwd(myloc1)
  QAstack <- stack("/Band_11.tif")
  #Create mask....
  QAstack[QAstack<128] <- NA
  print("mask created")
  #This is probably not right but anything > 128 is ok I think
  writeRaster(QAstack, "/QA_2015_2016.tif")
  B3stack <- stack("/Band_3.tif")
  B4stack <- stack("/Band_4.tif")
  B5stack <- stack("/Band_5.tif")
  B6stack <- stack("/Band_6.tif")
  B7stack <- stack("/Band_7.tif")
  B9stack <- stack("/Band_9.tif")
  B10stack <- stack("/Band_10.tif")
  print("masking raster stacks")
  B4 <- mask(B4stack, QAstack)
  B5 <- mask(B5stack, QAstack)
  B6 <- mask(B6stack, QAstack)
  B7 <- mask(B7stack, QAstack)
  B3 <- mask(B3stack, QAstack)
  B9 <- mask(B9stack, QAstack)
  B10 <- mask(B10stack, QAstack)
  print("writing masked rasters")
  writeRaster(B4, "/Band_4_Masked.tif")
  writeRaster(B5, "/Band_5_Masked.tif")
  writeRaster(B6, "/Band_6_Masked.tif")
  writeRaster(B7, "/Band_7_Masked.tif")
  writeRaster(B3, "/Band_3_Masked.tif")
  writeRaster(B9, "/Band_9_Masked.tif")
  writeRaster(B10, "/Band_10_Masked.tif")
  gc()
  return("Done")

}  
#Function to perform bandmath Operations
Bandmath <- function(tile){
  B3 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile, "_Band_3_Masked.tif", sep=""), overwrite=TRUE)
  B4 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/", tile, "_Band_4_Masked.tif", sep=""), overwrite=TRUE)
  B5 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/", tile, "_Band_5_Masked.tif", sep=""), overwrite=TRUE)
  B6 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/", tile, "_Band_6_Masked.tif", sep=""), overwrite=TRUE)
  B7 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/", tile, "_Band_7_Masked.tif", sep=""), overwrite=TRUE)
  
  NDVI <- ((B5 - B4)/(B5 + B4))
  writeRaster(NDVI, paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile,"_NDVI.tif", sep=""), overwrite=TRUE)
  STI <- (B6/B7)
  writeRaster(STI, paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile,"_STI.tif", sep=""), overwrite=TRUE)
  SINDRI <- ((B6 - B7)/(B6 + B7))
  writeRaster(SINDRI, paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile,"_SINDRI.tif", sep=""), overwrite=TRUE)
  #Important - before calculations....
  print("scalingNDVI")
  NDVI[NDVI>1] <- NA
  NDVI[NDVI<0] <- NA
  
  writeRaster(NDVI, paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile, "_NDVI_scaled.tif", sep=""), overwrite=TRUE)
  print("Bandmath Calcs")
  #For the band medians - just get the median of the first 8 observations------#changing from 2-10 (NDVI lag)
  fun8 <- function(x){median(x, na.rm=TRUE)}
  B3_med <- app(B3[[2:10]], fun8)
  B5_med <- app(B5[[2:10]], fun8)
  B6_med <- app(B6[[2:10]], fun8)
  NDVI_med <- app(NDVI[[2:10]], fun8)
  NDVI_mean <- app(NDVI[[2:10]], mean_na)
  NDVI_max <- app(NDVI[[2:10]], fun=function(x){max(x, na.rm=TRUE)})
  NDVI_min <- app(NDVI[[2:10]], fun=function(x){max(x, na.rm=TRUE)})
  NDVI_fullmax <- app(NDVI, fun=function(x){max(x, na.rm=TRUE)})
  NDVI_amp <- NDVI_fullmax-NDVI_max
  NDVI_ratio <- NDVI_med/NDVI_fullmax
  STI_med <- app(STI[[2:10]], fun8)
  SINDRI_med <- app(SINDRI[[2:10]], fun8)
  ### A much (> 100 times) faster approach is to directly use 
  ### linear algebra and pre-compute some constants
  
  which.max.na <- function(x, ...){
    max_idx <- which.max(x)
    ifelse(length(max_idx)==0,return(NA),return(max_idx))
    
  }
  print("GDD calc")
  GDD <-app(NDVI, which.max.na)
  
  
  inputs <- c(B3_med, B5_med, B6_med, NDVI_med, NDVI_mean, NDVI_max, NDVI_min, NDVI_fullmax, NDVI_amp, NDVI_ratio, GDD, SINDRI_med, STI_med)
  inputs
  print(str(inputs))
  names(inputs) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                     "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")
  writeRaster(inputs, paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/",tile, "_input_stack.tif"), overwrite=TRUE)
  gc()  
  print("done")
  
}
#list <- c("16SDH", "16SFH", "16SDJ", "16SEJ", "16SEH", "16SFJ", "16TFK", "16TEK", "16TDK", "16TDL", "16TEL", "16TFL")
Process_L30(x="16SDH", y="2015")
Process_L30(x="16SDH", y="2016")
write_bandstacks(x="16SDH", d="2016", c="2015")

Process_L30(x="16SEH", y="2015")
Process_L30(x="16SEH", y="2016")
write_bandstacks(x="16SEH", d="2016", c="2015")

Process_L30(x="16TDK", y="2015")
Process_L30(x="16TDK", y="2016")
write_bandstacks(x="16TDK", d="2016", c="2015")

Process_L30(x="16TEL", y="2015")
Process_L30(x="16TEL", y="2016")
write_bandstacks(x="16TEL", d="2016", c="2015")


#get average of first 5 observations (Oct Nov)

plot(NDVI, zlim=c(0,1))
plot(STI)
plot(SINDRI)
plot(B3stack)
plot(B4stack)
plot(B5stack)
plot(B6stack)
plot(B7stack)
plot(B9stack)
plot(B10stack)




rasterOptions(maxmemory = 1e+09, progress="text", overwrite=TRUE, chunksize=1e10)
list <- c("16SDH", "16SFH", "16SDJ", "16SEJ", "16SEH", "16SFJ", "16TFK", "16TEK", "16TDK", "16TDL", "16TEL", "16TFL")
#lapply(list[8:12], Bandmath)
Bandmath("16TDL")
Bandmath("16TEL")
Bandmath("16TEK")

#Format Windshield surveys by county----
format_windshield1 <-function(county){
  tst <- as.data.frame(shapefile(paste("/Volumes/G-RAID_Thunderbolt3/IN cover crop transect data/Buffer and transect data joins", county, paste(county, "Fall2015_Join.shp", sep="_"), sep="/")))
  tst <- tst[tst$R_Lon != 0,]
  #tst <- tst[tst$R_lon != 0,]
  print(head(tst))
  xy <- tst[,c(9,8)]
  print(head(xy))
  spdf <- SpatialPointsDataFrame(coords=xy, data=tst, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(as.list(c(xy, spdf)))
}
format_windshield4 <-function(county){
  tst <- as.data.frame(shapefile(paste("/Volumes/G-RAID_Thunderbolt3/IN cover crop transect data/Buffer and transect data joins", county, paste(county, "Fall2015_Join.shp", sep="_"), sep="/")))
  tst <- tst[tst$R_Lon != 0,]
  print(head(tst))
  #tst <- tst[tst$R_lon != 0,]
  xy <- tst[,c(10,9)]
  print(head(xy))
  spdf <- SpatialPointsDataFrame(coords=xy, data=tst, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(as.list(c(xy, spdf)))
}
format_windshield2 <-function(county){
  tst <- as.data.frame(shapefile(paste("/Volumes/G-RAID_Thunderbolt3/IN cover crop transect data/Buffer and transect data joins", county, paste(county, "Fall2015_Join.shp", sep="_"), sep="/")))
  #tst <- tst[tst$R_Lon != 0,]
  print(head(tst))
  tst <- tst[tst$R_lon != 0,]
  xy <- tst[,c(8,7)]
  print(head(xy))
  spdf <- SpatialPointsDataFrame(coords=xy, data=tst, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(as.list(c(xy, spdf)))
}
format_windshield3 <-function(county){
  tst <- as.data.frame(shapefile(paste("/Volumes/G-RAID_Thunderbolt3/IN cover crop transect data/Buffer and transect data joins", county, paste(county, "Fall2015_Join.shp", sep="_"), sep="/")))
  print(head(tst))
  #tst <- tst[tst$R_Lon != 0,]
  tst <- tst[tst$R_lon != 0,]
  xy <- tst[,c(9,8)]
  print(head(xy))
  spdf <- SpatialPointsDataFrame(coords=xy, data=tst, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(as.list(c(xy, spdf)))
}

Posey  <- format_windshield1("Posey")
Benton <- format_windshield2("Benton")
Gibson <- format_windshield3("Gibson")
Warren <- format_windshield1("Warren")
White <- format_windshield4("White")

tomerge <- as.list(list.files(path="/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands", recursive = TRUE, pattern = "*input_stack", full.names = TRUE))
str(tomerge)
tomerge <- tomerge[-3]

e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)
indiana <- raster(e)
proj4string(indiana) <- CRS('+init=epsg:4269')
writeRaster(indiana, file="IndianaInput.tif", format="GTiff")
mosaic_rasters(gdalfile=tomerge, dst_dataset="IndianaInput.tif", of="GTiff", output.vrt=NULL, separate=TRUE, verbose=TRUE, output_Raster = TRUE, CHECK_DISK_FREE_SPACE=FALSE)
gdalinfo("IndianaInput.tif")
#This file is giant-------
gc()
In <-terra::rast("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/IndianaInput.tif", NAvalue="-Inf")
#Mosaic didn't work quite like I thought it did. There's 169 layers here (13 x 13). So layers 1-13 are for 1 tile, 14-26 are for the second, etc. This means we're going to need to do some band math
#By layer
#In[In==-Inf] <- NA
#Spit apply combine
Set1 <- In[[1:13]]
Set2 <- In[[14:26]]
Set3 <- In[[27:39]]
Set4 <- In[[40:52]]
Set5 <- In[[53:65]]
Set6 <- In[[65:78]]
Set7 <- In[[79:91]]
Set8 <- In[[92:104]]
Set9 <- In[[105:117]]
Set10 <- In[[118:130]]
Set11 <- In[[131:143]]
Set12 <- In[[144:156]]

gc()
clamp(Set1, -1, 100)
Set1[Set1==-Inf] <- NA
coSet3[Set3==-Inf] <- NA
Set4[Set4==-Inf] <- NA
Set5[Set5==-Inf] <- NA
Set6[Set6==-Inf] <- NA
Set7[Set7==-Inf] <- NA
Set8[Set8==-Inf] <- NA
Set9[Set9==-Inf] <- NA
Set10[Set10==-Inf] <- NA
Set11[Set11==-Inf] <- NA
Set12[Set12==-Inf] <- NA

Set1[Set1==Inf] <- NA
function(x) {
  mean(x,na.rm=T)
}
plot(In[[2]])
min(In[[2]])

t16SDH<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDH_input_stack.tif")
names(t16SDH) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")
#16SCH<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SCDH_input_stack.tif")
#16SDG <- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDG_input_stack.tif")
t16TDK<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16TDK _input_stack.tif")
names(t16TDK) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")
t16TDL<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16TDL _input_stack.tif")
names(t16TDL) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")
t16TEL<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16TEL_input_stack.tif")
names(t16TEL) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")

Posey_Rf_input <- cbind(as.data.frame(Posey[[3]]), extract(t16SDH, (cbind(Posey[[1]] , Posey[[2]])), buffer=40, fun=median))
Gibson_Rf_input <- cbind(as.data.frame(Gibson[[3]]), extract(t16SDH, (cbind(Gibson[[1]] , Gibson[[2]])), buffer=40, fun=median))


Warren_coords <- cbind(Warren[[1]], Warren[[2]])
Warren_Rf_input <- cbind(as.data.frame(Warren[[3]]), extract(t16TDK, Warren_coords, buffer=40, fun=median))

White_coords <- cbind(White[[1]], White[[2]])
origin(t16TEL) <- 0
origin(t16TDL) <- 0
White_raster <- raster::mosaic(t16TEL, t16TDL, fun=median)
t16TEL
t16TDL
plot(t16TEL)
plot(t16TDL)
t16TEL
t16TDL
#White_ext <- rbind(extract(t16TDL, White_coords, buffer=40, fun=median), extract(t16TEL, White_coords, buffer=40, fun=median))
White_Rf_input <- cbind(as.data.frame(White[[3]]), extract(t16TEL, White_coords, buffer=10, fun=median))
#Benton kinda messed up - hold for validation?
Benton_Rf_input <- cbind(as.data.frame(Benton[[3]]), extract(t16TDK, (cbind(Benton[[1]] , Benton[[2]])), fun=median))
head(Posey_Rf_input)
Posey_Rf_input <- Posey_Rf_input[c("R_Lat", "R_Lon", "Prev_Crp", "Fall_Tilla", "Cover_Crop", "CC_Method", "B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")]

Posey_Rf_input$county <- "Posey"

head(Benton_Rf_input)
Benton_Rf_input <- Benton_Rf_input[c("R_lat", "R_lon", "Prev_Crp", "Fall_Tilla", "Cover_Crop", "CC_Method", "B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                                     "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")]
Benton_Rf_input <- rename(Benton_Rf_input, R_Lat=R_lat, R_Lon=R_lon)
Benton_Rf_input$county <- "Benton"

head(Gibson_Rf_input)
Gibson_Rf_input <- Gibson_Rf_input[c("R_lat", "R_lon", "Prev_Crp", "Fall_Tilla", "Cover_Crop", "CC_Method", "B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                                     "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")]

Gibson_Rf_input <- rename(Gibson_Rf_input, R_Lat=R_lat, R_Lon=R_lon)

Gibson_Rf_input$county <- "Gibson"

head(Warren_Rf_input)
Warren_Rf_input <- Warren_Rf_input[c("R_Lat", "R_Lon", "Prev_Crp", "Fall_Tilla", "Cover_Crop", "CC_Method", "B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                                     "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")]

Warren_Rf_input$county <- "Warren"

#rename(Warren_Rf_input, R_Lat=R_lat, R_Lon=R_lon)

head(White_Rf_input)
White_Rf_input <- White_Rf_input[c("R_Lat", "R_Lon", "Prev_Crp", "Fall_Tilla", "Cover_Crop", "CC_Method", "B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med")]
White_Rf_input$county <- "White"
#rename(White_Rf_input, R_Lat=R_lat, R_Lon=R_lon)

All_counties_input <- do.call("rbind", list(Posey_Rf_input, Benton_Rf_input, Gibson_Rf_input, Warren_Rf_input, White_Rf_input))
write.csv(All_counties_input, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/RS_input_all.csv")


#print(model)


