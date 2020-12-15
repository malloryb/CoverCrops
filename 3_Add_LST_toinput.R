#Landsat ST for Cover crops
#Helpful! https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/landsat-data-in-r-geotiff/
setwd("/Volumes/G-RAID_Thunderbolt3/Bulk Order Indiana Provisional Surface Temp/U.S. Landsat 4-8 ARD/")
crs <- "+proj=longlat +datum=WGS84 +no_defs +datum=WGS84 " # desired GeoTIFF SRS
fromCRS <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=IAU76 +units=m +no_defs"
tst <- raster("LC08_CU_021007_20151028_20181210_C01_V01_ST/LC08_CU_021007_20151028_20181210_C01_V01_ST.tif")
tst2 <- projectRaster(tst, crs=crs)
tst3 <- ((tst2/10)-273.15) #gotta rescale and convert to celcius
plot(tst3)

process_LandsatST <- function(x, tile=y){
  listQA <- list.files(path="/Volumes/G-RAID_Thunderbolt3/Bulk Order Indiana Provisional Surface Temp/U.S. Landsat 4-8 ARD/", recursive=TRUE, pattern="*STQA.tif", full.names=TRUE)
  tst <- raster::stack(x)
  crs <- "+proj=longlat +datum=WGS84 +no_defs +datum=WGS84 " # desired GeoTIFF SRS
  cloud <- listQA[grep(y, listQA)]
  print(cloud)
  cloud <- (raster::stack(cloud))
  print(ext(cloud))
  ext1 <- extent(tst)
  extent(cloud) <- ext1
  print("recoding QA raster")
  #This is for the level 1 QA, but not for our QA for the ST Product
  #[!cloud==2720 |  !cloud==2724 | !cloud==2728| !cloud==2732] <- NA
  cloud <- cloud*0.01 #scale factor: https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1330-LandsatSurfaceTemperature_ProductGuide-v2.pdf
  cloud[cloud > 7] <- NA
  #where y = the 6 digit landsat tile
  #read file
  print("masking raster")
  tst2 <- raster::mask(tst, cloud)
  print("projecting raster")
  tst2 <- projectRaster(tst2, crs=crs)
  print("projecting and converting to celsius")
  tst3 <- ((tst2/10)-273.15) #gotta rescale and convert to celcius
  #should mask out the clouds as well...
  return(tst3)
}

tst <- process_LandsatST(x=listST[grep("021010", listST)], tile="021010")
raster::plot(tst)
#raster::stack by tile and extract from there - figure out mosaicing later
listST <- list.files(path="/Volumes/G-RAID_Thunderbolt3/Bulk Order Indiana Provisional Surface Temp/U.S. Landsat 4-8 ARD", recursive = TRUE, pattern = "*_ST.tif", full.names = FALSE)
#Pull tile from filename
substr(listST, 9,14)
y <- "s021007"
s021007 <- raster::stack(lapply(listST[grep("021007", listST)], process_LandsatST, tile="021007"))
names(s021007)
writeRaster(s021007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021007.tif")
y <- "s021008"
s021008 <- raster::stack(lapply(listST[grep("021008", listST)], process_LandsatST))
writeRaster(s021008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021008.tif")
y <- "s021009"
s021009 <- raster::stack(lapply(listST[grep("021009", listST)], process_LandsatST))
writeRaster(s021009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021009.tif")
y <- "s021010"
s021010 <- raster::stack(lapply(listST[grep("021010", listST)], process_LandsatST, tile="021010"))
writeRaster(s021010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021010.tif")
y <- "s022007"
s022007 <- raster::stack(lapply(listST[grep("022007", listST)], process_LandsatST))
writeRaster(s022007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022007.tif")
y <- "s022008"
s022008 <- raster::stack(lapply(listST[grep("022008", listST)], process_LandsatST))
writeRaster(s022008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022008.tif")
y <- "s022009"
s022009 <- raster::stack(lapply(listST[grep("022009", listST)], process_LandsatST))
writeRaster(s022009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022009.tif")
y <- "s022010"
s022010 <- raster::stack(lapply(listST[grep("022010", listST)], process_LandsatST))
writeRaster(s022010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022010.tif")
y <- "s022011"
s022011 <- raster::stack(lapply(listST[grep("022011", listST)], process_LandsatST))
writeRaster(s022011, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022011.tif")
y <- "s023007"
s023007 <- raster::stack(lapply(listST[grep("023007", listST)], process_LandsatST))
writeRaster(s023007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023007.tif")
y <- "s023008"
s023008 <- raster::stack(lapply(listST[grep("023008", listST)], process_LandsatST))
writeRaster(s023008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023008.tif")
#s023008 <- raster::stack(lapply(listST[grep("023008", listST)], process_LandsatST))
y <- "s023009"
s023009 <- raster::stack(lapply(listST[grep("023009", listST)], process_LandsatST))
writeRaster(s023009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023009.tif")
y <- "s023010"
s023010 <- raster::stack(lapply(listST[grep("023010", listST)], process_LandsatST))
writeRaster(s023010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023010.tif")

SST_list[grep("021007", names(SST_list))]  
#For each tile, reproject and raster::stack up tifs (by date)

#extract #re-running - do we even need to scale? 
gc()
inputs <- read.csv("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/RS_input_all_12_8.csv")
inputs
coords.sp = as.data.frame(cbind(inputs$R_Lon, inputs$R_Lat))

s021007_fall <- s021007[[1:3,]]
#ts <- mean_na(scale(s021007_fall, center=TRUE, scale=TRUE))
ts <- mean(s021007_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s021007_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s021008_fall <- s021008[[1:3,]]
#ts <- mean_na(scale(s021008_fall, center=TRUE, scale=TRUE))
ts <- mean(s021008_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s021008_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s021009_fall <- s021009[[1:3,]]
#ts <- mean_na(scale(s021009_fall, center=TRUE, scale=TRUE))
ts <- mean(s021009_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s021009_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s021010_fall <- s021010[[1:3,]]
#ts <- mean_na(scale(s021010_fall, center=TRUE, scale=TRUE))
ts <- mean(s021010_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s021010_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s022007_fall <- s022007[[1:3,]]
#ts <- mean_na(scale(s022007_fall, center=TRUE, scale=TRUE))
ts <- mean(s022007_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s022007_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s022008_fall <- s022008[[1:3,]]
#ts <- mean_na(scale(s022008_fall, center=TRUE, scale=TRUE))
ts <- mean(s022008_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s022008_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s022009_fall <- s022009[[1:3,]]
#ts <- mean_na(scale(s022009_fall, center=TRUE, scale=TRUE))
ts <- mean(s022009_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s022009_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s022010_fall <- s022010[[1:3,]]
#ts <- mean_na(scale(s022010_fall, center=TRUE, scale=TRUE))
ts <- mean(s022010_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s022010_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s022011_fall <- s022011[[1:3,]]
#ts <- mean_na(scale(s022011_fall, center=TRUE, scale=TRUE))
ts <- mean(s022011_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s022011_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s023007_fall <- s023007[[1:3,]]
#ts <- mean_na(scale(s023007_fall, center=TRUE, scale=TRUE))
ts <- mean(s023007_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s023007_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s023008_fall <- s023008[[1:3,]]
#ts <- mean_na(scale(s023008_fall, center=TRUE, scale=TRUE))
ts <- mean(s023008_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s023008_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s023009_fall <- s023009[[1:3,]]
#ts <- mean_na(scale(s023009_fall, center=TRUE, scale=TRUE))
ts <- mean(s023009_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s023009_fall"=extract(ts, coords.sp, buffer=10, fun=median))

s023010_fall <- s023010[[1:3,]]
#ts <- mean_na(scale(s023010_fall, center=TRUE, scale=TRUE))
ts <- mean(s023010_fall, na.rm=TRUE)
inputs <- cbind(inputs, "s023010_fall"=extract(ts, coords.sp, buffer=10, fun=median))

head(inputs)
str(inputs)
inputs_LST <- inputs[,23:34]
str(inputs_LST)

#inputs$LST_scaled <- rowMeans(inputs_LST, na.rm=TRUE)

inputs$LST_unscaled <- rowMeans(inputs_LST, na.rm=TRUE)

inputs_write <- subset(inputs, select = -c(22:34)) %>%
  relocate(LST_unscaled, .before=county)


write.csv(inputs_write, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/inputs_with_LST_unscaled_12_8.csv")
