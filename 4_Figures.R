pacman::p_load(raster, ggplot2, gdalUtils, terra, dplyr, rgdal, sp, rhdf5, rlist, randomForest, caTools, caret, maptools, e1071)

#Figure A)
#Phenology of the different types of pixels
#Pull input raster for Posey

#Input NDVI
raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDH_NDVI_scaled.tif")
#Input LST 
raster::stack("/Volumes/G-RAID_Thunderbolt3/Bulk Order Indiana Provisional Surface Temp/U.S. Landsat 4-8 ARD/LC08_CU_023009_20160322_20181202_C0")
#Input SINDRI
raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDH_SINDRI.tif")
