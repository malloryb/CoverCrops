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
  return(tst3)
}

listST <- list.files(path="/Volumes/G-RAID_Thunderbolt3/Bulk Order Indiana Provisional Surface Temp/U.S. Landsat 4-8 ARD", recursive = TRUE, pattern = "*_ST.tif", full.names = FALSE)

s021007 <- process_LandsatST(x=listST[grep("021007", listST)], tile="021007")
writeRaster(s021007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021007.tif")

s021008 <- process_LandsatST(x=listST[grep("021008", listST)], tile="021008")
writeRaster(s021008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021008.tif")

s021009 <- process_LandsatST(x=listST[grep("021009", listST)], tile="021009")
writeRaster(s021009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021009.tif")

s021010 <- process_LandsatST(x=listST[grep("021010", listST)], tile="021010")
writeRaster(s021010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021010.tif")

s022007 <- process_LandsatST(x=listST[grep("022007", listST)], tile="022007")
writeRaster(s022007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022007.tif")

s022008 <- process_LandsatST(x=listST[grep("022008", listST)], tile="022008")
writeRaster(s022008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022008.tif")

s022009 <- process_LandsatST(x=listST[grep("022009", listST)], tile="022009")
writeRaster(s022009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022009.tif")

s022010 <- process_LandsatST(x=listST[grep("022010", listST)], tile="022010")
writeRaster(s022010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022010.tif")

s022011 <- process_LandsatST(x=listST[grep("022011", listST)], tile="022011")
writeRaster(s022011, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s022011.tif")

s023007 <- process_LandsatST(x=listST[grep("023007", listST)], tile="023007")
writeRaster(s023007, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023007.tif")

vs023008 <- process_LandsatST(x=listST[grep("023008", listST)], tile="023008")
writeRaster(s023008, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023008.tif")

s023009 <- process_LandsatST(x=listST[grep("023009", listST)], tile="023009")
writeRaster(s023009, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023009.tif")

s023010 <- process_LandsatST(x=listST[grep("023010", listST)], tile="023010")
writeRaster(s023010, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s023010.tif")

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
