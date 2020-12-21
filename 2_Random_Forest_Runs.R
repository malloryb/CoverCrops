#"Clean" Start to random forest model runs
pacman::p_load(raster, gdalUtils, terra, dplyr, rgdal, sp, rhdf5, rlist, randomForest, caTools, caret, maptools, e1071, caretEnsemble, lattice, gridExtra)

mean_na_x <- function(x) {
  e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)
  r2 <- raster(e)
  res(r2) <- c(0.00356, 0.00266)
  crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs"
  x <- stack(x)
  x <- x[[1:3]]
  y <- mean(x,na.rm=T)
  y <- extend(y, e)
  y <- resample(y,r2)
  print(res(y))
  return(y)
}

#First load all the stacks 
LSTmerged <- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/Input_LST.tif")
t16SDH<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16SDH _input_stack.tif")
w <- raster::crop(LSTmerged, t16SDH)
w <- raster::resample(w, t16SDH)
t16SDH <- raster::addLayer(t16SDH,w)
names(t16SDH) <-  c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                    "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax", "LST")

t16SEH<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16SEH _input_stack.tif")
x <- crop(LSTmerged, t16SEH)
x <- resample(x, t16SEH)
t16SEH <- raster::addLayer(t16SEH, x)
names(t16SEH) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax", "LST")

t16TDK<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TDK _input_stack.tif")
y <- crop(LSTmerged, t16TDK)
y <- resample(y, t16TDK)
t16TDK <- raster::addLayer(t16TDK, y)
names(t16TDK) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax", "LST")

t16TEL<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TEL _input_stack.tif")
z <- crop(LSTmerged, t16TEL)
z <- resample(z, t16TEL)
t16TEL <- raster::addLayer(t16TEL, z)
names(t16TEL) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax", "LST")

t16TDL<- raster::stack("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/ 16TDL _input_stack.tif")
a <- crop(LSTmerged, t16TDL)
a <- resample(a, t16TDL)
t16TDL <- raster::addLayer(t16TDL, a)
names(t16TDL) <- c("B3_med", "B5_med", "B6_med", "NDVI_med", "NDVI_mean", "NDVI_max", "NDVI_min", "NDVI_fullmax", "NDVI_amp", 
                   "NDVI_ratio", "GDD", "SINDRI_med", "STI_med", "B9_med", "B10_med", "therm_ratio", "B10_fullmax", "LST")

#input data csv------
gc()
All_counties_input <- read.csv("/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/input_all_12_11.csv")

#For today: 12/11, we're also going to add the proper LST data from the Greek online LST calculator
#Double Check 1 and 2 are the lat longs
head(All_counties_input)
All_counties_input$Cover_Crop <- as.factor(All_counties_input$Cover_Crop)
levels(All_counties_input$Cover_Crop)
plyr::count(All_counties_input$Cover_Crop)

levels(as.factor(All_counties_input$Fall_Tilla))
plyr::count(All_counties_input$Fall_Tilla)

levels(as.factor(All_counties_input$CC_Method))
plyr::count(All_counties_input$CC_Method)

#Two class model-------
#Cover Crop Presence Absence Category: 0 = no cover crop,  2 = "Cover Crop" 
#Cover Crop Specific Category: 0 = nothing (presumably conventional tillage), 1= 'No-Till', 2 = "Cover crop".  

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

counts <- (subset(All_counties_input, county=="Warren"))
plyr::count(counts$Cover_Crop_PA)
#Posey: 52/132
#Gibson: 55/188
#White: 5/190
#Benton:6/154
#Warren: 14/193
#random forest model develop
#Gonna try github method from here: https://gist.github.com/hakimabdi/720f1481af9eca0b7b97d9856052e0e2
# Split the data frame into 70-30 by class
#First grouping------
# Create model weights (they sum to one)
set.seed(400)
str(All_counties_input)
All_counties_input <- upSample(All_counties_input, All_counties_input$Cover_Crop_PA)
sample <-sample.split(All_counties_input$Cover_Crop_PA, SplitRatio = 0.8)
trn=subset(All_counties_input, sample==TRUE)
eva=subset(All_counties_input, sample==FALSE)
eva.sp = SpatialPointsDataFrame(coords = cbind(eva$R_Lon, eva$R_Lat), data = eva, proj4string = crs("+proj=longlat +datum=WGS84 +no_defs"))
# Set up a resampling method in the model training process
tc <- trainControl(method = "rf", # repeated cross-validation of the training data
                     number = 10, # number of folds
                   repeats = 10, # number of repeats
                   allowParallel = FALSE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE) # view the training iterations

# Generate a grid search of candidate hyper-parameter values for inclusion into the models training
# These hyper-parameter values are examples. You will need a more complex tuning process to achieve high accuracies
# For example, you can play around with the parameters to see which combinations gives you the highest accuracy. 
rf.grid <- expand.grid(mtry=1:19) # number of variables available for splitting at each tree node

## Begin training the models. It took my laptop 8 minutes to train all three algorithms
# Train the random forest model
head(trn)
#THIS IS WHERE WE CHANGE ROWS 
rf_modelLST <- caret::train(x = trn[,(9:26)], y = as.factor(trn$Cover_Crop_PA),
                         method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelSWIR <- caret::train(x = trn[,(9:21)], y = as.factor(trn$Cover_Crop_PA),
                         method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelVISNir <- caret::train(x = trn[,(9:19)], y = as.factor(trn$Cover_Crop_PA),
                         method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelNDVI <- caret::train(x = trn[12], y = as.factor(trn$Cover_Crop_PA),
                         method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)

rf_modelLST
rf_modelSWIR
rf_modelVISNir
rf_modelNDVI

varImp(rf_modelNDVI)
varImp(rf_modelVISNir)
varImp(rf_modelSWIR)
varImp(rf_modelLST)
save(rf_modelLST,file = "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/RandomForest_LST_PA_12_16.RData")

#save(rf_model,file = "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/RandomForest_LST_PA.RData")
#random forest model evalaute
## Apply the models to data. It took my imac 16 minutes to apply the random forest model
# Apply the random forest model to the HLS data-----
gc()
rf_prediction1 = raster::predict(t16SDH, model=rf_modelLST)
rf_prediction1_2 = raster::predict(t16SDH, model=rf_modelSWIR)
rf_prediction1_3 = raster::predict(t16SDH, model=rf_modelVISNir)
rf_prediction1_4 = raster::predict(t16SDH, model=rf_modelNDVI)

writeRaster(rf_prediction1, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/LST_model_PA_12_15.tif")

raster::plot(rf_prediction1)
library(ggplot2)
library(reshape2)
raster::hist(rf_prediction1)
rasterVis::histogram(rf_prediction1)
rasterVis::histogram(rf_prediction2)
rasterVis::histogram(rf_prediction3)
rasterVis::histogram(rf_prediction4)

varImpPlot(rf_modelLST)
rf_prediction2= raster::predict(t16TDK, model=rf_modelLST)
rf_prediction3 = raster::predict(t16TDL, model=rf_modelLST)
rf_prediction4 = raster::predict(t16TEL, model=rf_modelLST)

varImp(rf_modelLST)
# Custom function ---------------------------------------------------------
# The function requires list of models as input and is used in for loop 
plot_importance <- function(importance_list, imp, algo_names) {
  importance <- importance_list[[imp]]$importance
  model_title <- algo_names[[imp]]
  if (ncol(importance) < 2) { # Plot dotplot if dim is ncol < 2
    importance %>%
      as.matrix() %>%
      dotplot(main = model_title)
  } else { # Plot heatmap if ncol > 2
    importance %>%
      as.matrix() %>%
      levelplot(xlab = NULL, ylab = NULL, main = model_title, scales = list(x = list(rot = 45)))
  }
}
my_list_of_models <- c(rf_modelLST, rf_modelLST2)
importance <- lapply(my_list_of_models, varImp)

# Plotting variable immportance -------------------------------------------
# Create second loop to go over extracted importance and plot it using plot()
importance_plots <- list()
for (imp in seq_along(importance)) {
  # importance_plots[[imp]] <- plot(importance[[imp]])
  importance_plots[[imp]] <- plot_importance(importance_list = importance, imp = imp, algo_names = names(my_list_of_models))
}

importance_plots
# Multiple plots at once
do.call("grid.arrange", c(importance_plots))


rf1In <- extend(rf_prediction1, e)  
rf2In <- extend(rf_prediction2, e)  
rf3In <- extend(rf_prediction3, e)  
rf4In <- extend(rf_prediction4, e)  

overlay(rf1In, rf2In, rf3In, rf4In, fun=mean_na)
e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)

raster::plot(rf_prediction1)
#plot(rf_prediction2)
#plot(rf_prediction3)
#plot(rf_prediction4)

#rf_prediction_prob = raster::predict(t16SDH, model=rf_model, type="prob")

# Convert the evaluation data into a spatial object using the X and Y coordinates and extract predicted values
# Create an error matrix for each of the classifiers
rf_Eval1 = extract(rf_prediction1, eva.sp)
rf_Eval1_2 = extract(rf_prediction1_2, eva.sp)
rf_Eval1_3 = extract(rf_prediction1_3, eva.sp)
rf_Eval1_4 = extract(rf_prediction1_4, eva.sp)

rf_Eval2 = extract(rf_prediction2, eva.sp)
rf_Eval3 = extract(rf_prediction3, eva.sp)
rf_Eval4 = extract(rf_prediction4, eva.sp)
#e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)
#r2 <- raster(e)
#res(r2) <- c(0.00356, 0.00266)
#crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs"
#f <- extend(rf_prediction1, e)
#f <- resample(f, r2)
#d <- extend(rf_prediction2, e)
#d <- resample(d,r2)
#c <- extend(rf_prediction3, e)
#c <- resample(c,r2)
#b <- extend(rf_prediction4, e)
#b <- resample(b,r2)
#plot(calc(stack(f,d,c,b), mean_na))

rf_errorM1 = confusionMatrix(recode(as.factor(rf_Eval1), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA),positive="ClassC")
rf_errorM1_2 = confusionMatrix(recode(as.factor(rf_Eval1_2), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM1_3 = confusionMatrix(recode(as.factor(rf_Eval1_3), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM1_4 = confusionMatrix(recode(as.factor(rf_Eval1_4), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM1
rf_errorM1_2
rf_errorM1_3
rf_errorM1_4


# From https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package/42940553
draw_confusion_matrix <- function(cm) {
  
  total <- sum(cm$table)
  res <- as.numeric(cm$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
} 
draw_confusion_matrix(rf_errorM1)
draw_confusion_matrix(rf_errorM1_2)
draw_confusion_matrix(rf_errorM1_3)
draw_confusion_matrix(rf_errorM1_4)

rf_errorM2 = confusionMatrix(recode(as.factor(rf_Eval2), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM3 = confusionMatrix(recode(as.factor(rf_Eval3), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM4 = confusionMatrix(recode(as.factor(rf_Eval4), "1"= "Class0", "2"="ClassC"),as.factor(eva$Cover_Crop_PA), positive="ClassC")
rf_errorM1
rf_errorM2
rf_errorM3
rf_errorM4

gc()
#Second grouping- 3 class model-----
levels(All_counties_input$Cover_Crop_Specific)
All_counties_input <- upSample(All_counties_input, All_counties_input$Cover_Crop_Specific)
sample <-sample.split(All_counties_input$Cover_Crop_Specific, SplitRatio = 0.80)
trn=subset(All_counties_input, sample==TRUE)
eva=subset(All_counties_input, sample==FALSE)
levels(All_counties_input$Cover_Crop_Specific)
# Set up a resampling method in the model training process
tc <- trainControl(method = "repeatedcv", # repeated cross-validation of the training data
                   number = 10, # number of folds
                   repeats = 5, # number of repeats
                   allowParallel = FALSE, # allow use of multiple cores if specified in training
                   verboseIter = TRUE) # view the training iterations

# Generate a grid search of candidate hyper-parameter values for inclusion into the models training
# These hyper-parameter values are examples. You will need a more complex tuning process to achieve high accuracies
# For example, you can play around with the parameters to see which combinations gives you the highest accuracy. 
rf.grid <- expand.grid(mtry=1:10) # number of variables available for splitting at each tree node

## Begin training the models. It took my laptop 8 minutes to train all three algorithms
# Train the random forest model
#rf_model2 <- caret::train(x = trn[,(9:22)], y = as.factor(trn$Cover_Crop_Specific),
#                         method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)

rf_modelLST2 <- caret::train(x = trn[,(9:26)], y = as.factor(trn$Cover_Crop_Specific),
                            method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelSWIR2 <- caret::train(x = trn[,(9:21)], y = as.factor(trn$Cover_Crop_Specific),
                             method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelVISNir2 <- caret::train(x = trn[,(9:19)], y = as.factor(trn$Cover_Crop_Specific),
                               method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)
rf_modelNDVI2 <- caret::train(x = trn[12], y = as.factor(trn$Cover_Crop_Specific),
                             method = "rf", metric="Kappa", trainControl = tc, tuneGrid = rf.grid)

rf_modelLST2
rf_modelSWIR2
rf_modelNDVI2
rf_modelVISNir2


rf_model2
varImp(rf_modelLST)
save(rf_modelLST2,file = "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/2015_2016_Input_Bands/RandomForest_LST_Cat_12_16.RData")

#random forest model evalaute
## Apply the models to data. It took my imac 16 minutes to apply the random forest model
# Apply the random forest model to the HLS data
rf2_prediction1 = raster::predict(t16SDH, model=rf_modelLST2)
rf2_prediction1_2 = raster::predict(t16SDH, model=rf_modelSWIR2)
rf2_prediction1_3 = raster::predict(t16SDH, model=rf_modelVISNir2)
rf2_prediction1_4 = raster::predict(t16SDH, model=rf_modelNDVI2)

writeRaster(rf2_prediction1, "/Volumes/G-RAID_Thunderbolt3/HLS30_Indiana/LST_model_3class_12_11.tif")

rf2_prediction2= raster::predict(t16TDK, model=rf_model2)
rf2_prediction3 = raster::predict(t16TDL, model=rf_model2)
rf2_prediction4 = raster::predict(t16TEL, model=rf_model2)

#This is to plot everything on the same map (I'll want to fix up the rest of these too)
e <- extent(-88.09776,-84.784579,	37.771742, 41.760592)
r2 <- raster(e)
res(r2) <- c(0.00356, 0.00266)
crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs"
x <- extend(rf2_prediction1, e)
x <- resample(x, r2)
y <- extend(rf2_prediction2, e)
y <- resample(y,r2)
z <- extend(rf2_prediction3, e)
z <- resample(z,r2)
a <- extend(rf2_prediction4, e)
a <- resample(a,r2)

tst <- stack(x,y,z,a)
plot(calc(tst,mean_na))

#rf_prediction_prob = raster::predict(t16SDH, model=rf_model, type="prob")

# Convert the evaluation data into a spatial object using the X and Y coordinates and extract predicted values
# Create an error matrix for each of the classifiers
length(eva$Cover_Crop_Specific)
eva.sp = SpatialPointsDataFrame(coords = cbind(eva$R_Lon, eva$R_Lat), data = eva, proj4string = crs("+proj=longlat +datum=WGS84 +no_defs"))
length(eva.sp)
rf2_Eval1 = extract(rf2_prediction1, eva.sp)
rf2_Eval1_2 = extract(rf2_prediction1_2, eva.sp)
rf2_Eval1_3 = extract(rf2_prediction1_3, eva.sp)
rf2_Eval1_4 = extract(rf2_prediction1_4, eva.sp)

rf2_Eval2 = extract(rf2_prediction2, eva.sp)
rf2_Eval3 = extract(rf2_prediction3, eva.sp)
rf2_Eval4 = extract(rf2_prediction4, eva.sp)


rf2_errorM1 = confusionMatrix(recode(as.factor(rf2_Eval1), "1"= "Class0", "2"="ClassC", "3"="ClassT"),as.factor(eva$Cover_Crop_Specific))
rf2_errorM1_2 = confusionMatrix(recode(as.factor(rf2_Eval1_2), "1"= "Class0", "2"="ClassC", "3"="ClassT"),as.factor(eva$Cover_Crop_Specific))
rf2_errorM1_3 = confusionMatrix(recode(as.factor(rf2_Eval1_3), "1"= "Class0", "2"="ClassC", "3"="ClassT"),as.factor(eva$Cover_Crop_Specific))
rf2_errorM1_4 = confusionMatrix(recode(as.factor(rf2_Eval1_4), "1"= "Class0", "2"="ClassC", "3"="ClassT"),as.factor(eva$Cover_Crop_Specific))
draw_confusion_matrix(rf2_errorM1)
draw_confusion_matrix(rf2_errorM1_2)
draw_confusion_matrix(rf2_errorM1_3)
draw_confusion_matrix(rf2_errorM1_4)

rf2_errorM2 = confusionMatrix(as.factor(rf2_Eval2),as.factor(eva$Cover_Crop_Specific))
rf2_errorM3 = confusionMatrix(as.factor(rf2_Eval3),as.factor(eva$Cover_Crop_Specific))
rf2_errorM4 = confusionMatrix(as.factor(rf2_Eval4),as.factor(eva$Cover_Crop_Specific))
rf2_errorM1
rf2_errorM2
rf2_errorM3
rf2_errorM4
#rf_prediction_prob = raster::predict(t16SDH, model=rf_model, type="prob")


# Convert the evaluation data into a spatial object using the X and Y coordinates and extract predicted values
plot(rf_prediction2)
writeRaster(rf_prediction2,"/Volumes/Macintosh HD/Users/malbarn/Documents/CoverCrop_Project/rf_prediction_posey_cat2.tif", overwrite=TRUE) 
#writeRaster(rf_prediction_prob, "/Volumes/Macintosh HD/Users/malbarn/Documents/CoverCrop_Project/rf_prediction_prob_posey_cat2.tif")
plot(rf_prediction_prob)
varImp(rf_model2)
## Superimpose evaluation points on the predicted classification and extract the values
# random forest
rf_Eval = extract(rf_prediction2, eva.sp)

# Create an error matrix for each of the classifiers
rf_errorM2 = confusionMatrix(as.factor(rf_Eval),as.factor(eva$Cover_Crop_Specific))
rf_errorM2
#plots