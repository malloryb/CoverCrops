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
masking <- function(x){
  myloc1=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "2015-2016wy",  sep="/")
  print(myloc1)
  setwd(myloc1)
  myloc2=paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana", x, "QAstack.tif", sep="/")
  QAstack <- raster::stack("Band_11.tif")
  #Create mask....
  #This is probably not right but anything > 128 is ok I think
  QAstack[QAstack<128] <- NA
  print("mask created")
  writeRaster(QAstack, myloc2, overwrite=TRUE)
  B3stack <- raster::stack("Band_3.tif")
  B4stack <- raster::stack("Band_4.tif")
  B5stack <- raster::stack("Band_5.tif")
  B6stack <- raster::stack("Band_6.tif")
  B7stack <- raster::stack("Band_7.tif")
  B9stack <- raster::stack("Band_9.tif")
  B10stack <- raster::stack("Band_10.tif")
  print("masking raster stacks and writing masked rasters")
  B3 <- mask(B3stack, QAstack)
  writeRaster(B3, "Band_3_Masked.tif")
  gc()
  B4 <- mask(B4stack, QAstack)
  writeRaster(B4, "Band_4_Masked.tif")
  gc()
  B5 <- mask(B5stack, QAstack)
  writeRaster(B5, "Band_5_Masked.tif")
  gc()
  B6 <- mask(B6stack, QAstack)
  writeRaster(B6, "Band_6_Masked.tif")
  gc()
  B7 <- mask(B7stack, QAstack)
  writeRaster(B7, "Band_7_Masked.tif")
  gc()
  B9 <- mask(B9stack, QAstack)
  writeRaster(B9, "Band_9_Masked.tif")
  gc()
  B10 <- mask(B10stack, QAstack)
  gc()
  writeRaster(B10, "Band_10_Masked.tif")
  gc()
  return("Done")
}  
#Function to perform bandmath Operations
Bandmath <- function(tile){
  B3 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/",tile, "/2015-2016wy/", "Band_3_Masked.tif", sep=""), overwrite=TRUE)
  B4 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile, "/2015-2016wy/", "Band_4_Masked.tif", sep=""), overwrite=TRUE)
  B5 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile, "/2015-2016wy/", "Band_5_Masked.tif", sep=""), overwrite=TRUE)
  B6 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile,"/2015-2016wy/",  "Band_6_Masked.tif", sep=""), overwrite=TRUE)
  B7 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile,"/2015-2016wy/",  "Band_7_Masked.tif", sep=""), overwrite=TRUE)
  B9 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile, "/2015-2016wy/", "Band_9_Masked.tif", sep=""), overwrite=TRUE)
  B10 <- terra::rast(paste("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/", tile,"/2015-2016wy/",  "Band_10_Masked.tif", sep=""), overwrite=TRUE)
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
  gc()
  #For the band medians - just get the median of the first 8 observations------#changing from 2-10 (NDVI lag)
  fun8 <- function(x){median(x, na.rm=TRUE)}
  B3_med <- app(B3[[2:10]], fun8)
  B5_med <- app(B5[[2:10]], fun8)
  B6_med <- app(B6[[2:10]], fun8)
  B7_med <- app(B7[[2:10]], fun8)
  B9_med <- app(B9[[2:10]], fun8)
  B10_med <- app(B10[[2:10]], fun8)
  NDVI_med <- app(NDVI[[2:10]], fun8)
  print("all band medians finished")
  NDVI_mean <- app(NDVI[[2:10]], mean_na)
  NDVI_max <- app(NDVI[[2:10]], fun=function(x){max(x, na.rm=TRUE)})
  NDVI_min <- app(NDVI[[2:10]], fun=function(x){max(x, na.rm=TRUE)})
  NDVI_fullmax <- app(NDVI, fun=function(x){max(x, na.rm=TRUE)})
  print("calculating ratios")
  B10_fullmax <-app(B10, fun=function(x){max(x, na.rm=TRUE)})
  NDVI_amp <- NDVI_fullmax-NDVI_max
  NDVI_ratio <- NDVI_med/NDVI_fullmax
  therm_ratio <- (B10_med/B10_fullmax)
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
  
  
  inputs <- c(B3_med, B5_med, B6_med, NDVI_med, NDVI_mean, NDVI_max, NDVI_min, NDVI_fullmax, NDVI_amp, NDVI_ratio, GDD, SINDRI_med, STI_med, B9_med, B10_med, therm_ratio, B10_fullmax)
  inputs
  print(str(inputs))
  names(inputs) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                     "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax")
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

#masking
masking(x="16SDH")
masking(x="16SEH")
masking(x="16TDK")
masking(x="16TEL")



rasterOptions(maxmemory = 1e+09, progress="text", overwrite=TRUE, chunksize=1e10)
list <- c("16SDH", "16SFH", "16SDJ", "16SEJ", "16SEH", "16SFJ", "16TFK", "16TEK", "16TDK", "16TDL", "16TEL", "16TFL")
#lapply(list[8:12], Bandmath)
Bandmath("16SDH")
Bandmath("16SEH")
Bandmath("16TDK")
Bandmath("16TEL")

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

#Merge stuff - save for later
#tomerge <- as.list(list.files(path="/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands", recursive = TRUE, pattern = "*input_stack", full.names = TRUE))
#tomerge <- tomerge[-3]
#e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)
#indiana <- raster(e)
#proj4string(indiana) <- CRS('+init=epsg:4269')
#writeRaster(indiana, file="IndianaInput.tif", format="GTiff")
#mosaic_rasters(gdalfile=tomerge, dst_dataset="IndianaInput.tif", of="GTiff", output.vrt=NULL, separate=TRUE, verbose=TRUE, output_Raster = TRUE, CHECK_DISK_FREE_SPACE=FALSE)
#gdalinfo("IndianaInput.tif")

t16SDH<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16SDH _input_stack.tif")
names(t16SDH) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax")

#16SCH<- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SCDH_input_stack.tif")
#16SDG <- stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDG_input_stack.tif")
t16TDK<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TDK _input_stack.tif")
names(t16TDK) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax")

#t16TDL<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TDL _input_stack.tif")
#names(t16TDL) <-c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                  "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax")

t16TEL<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TEL _input_stack.tif")
names(t16TEL) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax")

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


