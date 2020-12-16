pacman::p_load(raster, ggplot2, gdalUtils, terra, dplyr, rgdal, sp, rhdf5, rlist, randomForest, caTools, caret, maptools, e1071)

#Figure A)
#Phenology of the different types of pixels
#Pull input raster for Posey

#Input NDVI
NDVIstack <- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDH_NDVI_scaled.tif")

#Input LST 
LSTstack <- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/s021010.tif")
#Input SINDRI
SINDRIstack <- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/16SDH_SINDRI.tif")
#Input data
All_counties_input <- read.csv("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/input_all_12_11.csv")
All_counties_input$Cover_Crop_Specific <- "Class0"
All_counties_input$Cover_Crop_Specific[All_counties_input$Fall_Tilla == "N" | All_counties_input$Fall_Tilla == "n"] <- "ClassT"
All_counties_input$Cover_Crop_Specific[All_counties_input$Cover_Crop != "N"] <- "ClassC"
All_counties_input$Cover_Crop_Specific <- as.factor(All_counties_input$Cover_Crop_Specific)
levels(as.factor(All_counties_input$Cover_Crop_Specific))

#recode cover crop: 
All_counties_input$Cover_Crop_PA <- "Class0"
All_counties_input$Cover_Crop_PA[All_counties_input$Cover_Crop == "A" | All_counties_input$Cover_Crop=="B" | All_counties_input$Cover_Crop=="BC"| All_counties_input$Cover_Crop=="C"| All_counties_input$Cover_Crop=="W" | All_counties_input$Cover_Crop=="G"] <- "ClassC"
All_counties_input$Cover_Crop_PA <- as.factor(All_counties_input$Cover_Crop_PA)
head(All_counties_input)
#All_counties_input <- subset(All_counties_input, select=-c(B9_med))
All_counties_input[sapply(All_counties_input, is.infinite)] <- NA
All_counties_input <- na.omit(All_counties_input)
plyr::count(All_counties_input$Cover_Crop_PA)
plyr::count(All_counties_input$Cover_Crop_Specific)


#Now, extract time series from the raster stacks
Cover <- subset(All_counties_input, county=="Posey" & Cover_Crop_Specific=="ClassC")
Conventional <- subset(All_counties_input, county=="Posey" & Cover_Crop_Specific=="Class0")
Residue <- subset(All_counties_input, county=="Posey" & Cover_Crop_Specific=="ClassT")
head(Cover)
head(Conventional)
head(Residue)
#Now, cbind everything together?
#Gonna scale within the time series of the year, right? 
#Try for one point first.

#For cover crops
coords.sptst = as.data.frame(cbind(Cover$R_Lon, Cover$R_Lat))
t1 <- as.data.frame(extract(NDVIstack, coords.sptst))
t2 <- as.data.frame(extract(LSTstack, coords.sptst))
t3 <- as.data.frame(extract(SINDRIstack, coords.sptst))

NDVIxCCseason <- (colMeans(t1, na.rm=TRUE))
NDVIxCCseason[is.nan(NDVIxCCseason)] <- NA


SWIRxCCseason <- (colMeans(t3, na.rm=TRUE))
SWIRxCCseason[is.nan(SWIRxCCseason)] <- NA


LSTxCCseason <- (colMeans(t2, na.rm=TRUE))
LSTxCCseason[is.nan(LSTxCCseason)] <- NA

#For crop residues/no-till 
coords.sptst2 = as.data.frame(cbind(Residue$R_Lon, Residue$R_Lat))
t4 <- as.data.frame(extract(NDVIstack, coords.sptst2))
t5 <- as.data.frame(extract(LSTstack, coords.sptst2))
t6 <- as.data.frame(extract(SINDRIstack, coords.sptst2))

NDVIxRseason <- (colMeans(t4, na.rm=TRUE))
NDVIxRseason[is.nan(NDVIxRseason)] <- NA


SWIRxRseason <- (colMeans(t6, na.rm=TRUE))
SWIRxRseason[is.nan(SWIRxRseason)] <- NA


LSTxRseason <- (colMeans(t5, na.rm=TRUE))
LSTxRseason[is.nan(LSTxRseason)] <- NA
#For conventional tillage
coords.sptst3 = as.data.frame(cbind(Conventional$R_Lon, Conventional$R_Lat))
t7 <- as.data.frame(extract(NDVIstack, coords.sptst3))
t8 <- as.data.frame(extract(LSTstack, coords.sptst3))
t9 <- as.data.frame(extract(SINDRIstack, coords.sptst3))

NDVIxCseason <- (colMeans(t7, na.rm=TRUE))
NDVIxCseason[is.nan(NDVIxCseason)] <- NA

SWIRxCseason <- (colMeans(t9, na.rm=TRUE))
SWIRxCseason[is.nan(SWIRxCseason)] <- NA

LSTxCseason <- (colMeans(t8, na.rm=TRUE))
LSTxCseason[is.nan(LSTxCseason)] <- NA


#Merge t1, t4, and t7 together....x should be 1-45, and y should be "Conventional", "Reisdue", and "Cover Crop

df <- data.frame(date=numeric(45), NDVICC=numeric(45), NDVIR=numeric(45), NDVICon=numeric(45))
df$date <- c(1:45)
df$NDVICC <- NDVIxCCseason
df$NDVIR <-NDVIxRseason
df$NDVICon <- NDVIxCseason

df2 <- data.frame(date=numeric(45), SWIRCC=numeric(45), SWIRR=numeric(45), SWIRCon=numeric(45))
df2$date <- c(1:45)
df2$SWIRCC <- SWIRxCCseason
df2$SWIRR <-SWIRxRseason
df2$SWIRCon <- SWIRxCseason

df3 <- data.frame(date=numeric(22), LSTCC=numeric(22), LSTR=numeric(22), LSTCon=numeric(22))
df3$date <- c(1:22)
df3$LSTCC <- LSTxCCseason
df3$LSTR <-LSTxRseason
df3$LSTCon <- LSTxCseason

NDVidata <- melt(df, id=c("date"))
ggplot(NDVidata[!is.na(NDVidata$value),], aes(x=date, y=value, color=variable, group=variable))+
  #geom_line()+
  geom_smooth(span=0.1)+
  theme_minimal()


SWIRdata <- melt(df2, id=c("date"))
ggplot(SWIRdata[!is.na(SWIRdata$value),], aes(x=date, y=value, color=variable, group=variable))+
  #geom_line()+
  geom_smooth(span=0.1)+
  theme_minimal()


LSTdata <- melt(df3, id=c("date"))
ggplot(LSTdata[!is.na(LSTdata$value),], aes(x=date, y=value, color=variable, group=variable))+
  geom_line()+
  #geom_smooth(span=0.3)+
  theme_minimal()

#Going to try setting conventional agriculture to zero here. 

df$NDVIR_diff <- df$NDVIR - df$NDVICon 
df$NDVICC_diff <- df$NDVICC - df$NDVICon 
NDVidata2 <- melt(df[,c("date", "NDVIR_diff", "NDVICC_diff")], id=c("date"))
ggplot(NDVidata2[!is.na(NDVidata2$value),], aes(x=date, y=value, color=variable, group=variable))+
  #geom_line()+
  ylab("NDVI")+
  geom_smooth(span=0.5, se=FALSE)+
  geom_hline(aes(yintercept=0))+
  theme_minimal(base_size=18)


df2$SWIRR_diff <- df2$SWIRR - df2$SWIRCon 
df2$SWIRCC_diff <- df2$SWIRCC - df2$SWIRCon 
SWIRdata2 <- melt(df2[,c("date", "SWIRR_diff", "SWIRCC_diff")], id=c("date"))
ggplot(SWIRdata2[!is.na(NDVidata2$value),], aes(x=date, y=value, color=variable, group=variable))+
  #geom_line()+
  ylab("SINDRI")+
  geom_smooth(span=0.5, se=FALSE)+
  geom_hline(aes(yintercept=0))+
  theme_minimal(base_size=18)


df3$LSTR_diff <- df3$LSTR - df3$LSTCon 
df3$LSTCC_diff <- df3$LSTCC - df3$LSTCon 
LSTdata2 <- melt(df3[,c("date", "LSTR_diff", "LSTCC_diff")], id=c("date"))
ggplot(LSTdata2[!is.na(NDVidata2$value),], aes(x=date, y=value, color=variable, group=variable))+
  #geom_line()+
  ylab("LST") +
  geom_smooth(span=0.5, se=FALSE)+
  geom_hline(aes(yintercept=0))+
  theme_minimal(base_size=18)



LST_Appears <- read.csv("/Volumes/G-RAID_Thunderbolt3/CC-more-MOD11A1-006-results.csv")
NDVI_Appears <- read.csv("/Volumes/G-RAID_Thunderbolt3/CC-more-MOD13Q1-006-results.csv")

LST_Appears$MOD11A1_006_LST_Day_1km
LST_Appears$Date
LST_Appears$Category

NDVI_Appears$Date
NDVI_Appears$MOD13Q1_006__250m_16_days_NDVI
NDVI_Appears$Category

ggplot(NDVI_Appears, aes(x=as.Date(Date), y=MOD13Q1_006__250m_16_days_NDVI, group=Category, color=Category))+
  geom_smooth(se=FALSE, span=0.2)+
  ylim(0,1)+
  xlab("Date")+
  ylab("NDVI")+
  theme_minimal(base_size=20)

ggplot(LST_Appears, aes(x=as.Date(Date), y=MOD11A1_006_LST_Day_1km, group=Category, color=Category))+
  geom_smooth(se=FALSE, span=0.1)+
  xlab("Date")+
  ylab("LST")+
  theme_minimal(base_size=20)