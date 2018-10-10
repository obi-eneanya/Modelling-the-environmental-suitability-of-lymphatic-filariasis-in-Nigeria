
## ~~~~~~~~~~~ Modelling occurrence (environmental suitability) of Lymphatic Filariasis in Nigeria ~~~~~~~~~~~~~~~ ##


# Set up a function to install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# List of packages
packages <- c("sp","raster","dismo","maptools","rgdal","proj4","ggplot2","gridExtra","cowplot","corrplot","biomod2","rJava","pROC","matrixStats","usdm","kernlab","ks","sm","gbm","mgcv","nlme","Metrics","tidyr")
ipak(packages)

# Import shapefiles
shapefiles <- setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Shapefiles")
Nigeria.ADM0 <- subset(Africa_ADM0, Africa_ADM0$ADM0_NAME == "Nigeria")
Nigeria.ADM1 <- subset(Africa_ADM1, Africa_ADM1$ADM0_NAME == "Nigeria")

plot(Nigeria.ADM0)
plot(Nigeria.ADM1, add= TRUE)

#load LF survey data#
NG.LF <- read.csv("NG.LF.csv") #this does not include the covariates

# Set up occurrence when at least 1 LF case is diagnosed by any diagnosis method.
NG.LF$Occurrence[NG.LF$Occurrence > 0] <- "1"
NG.LF$Occurrence[NG.LF$Occurrence == 0] <- "0"

#changing "1" and "0" from character to numeric
NG.LF <- transform(NG.LF, Occurrence = as.numeric(Occurrence)) 
class(NG.LF$Occurrence)
barplot(table(NG.LF$Occurrence)) # 324 absences and 1054 presences #old value
names(NG.LF$Occurrence)
summary(NG.LF)

#Exclude all post-MDA data
NG.LF <- subset(NG.LF, NG.LF$PoT == c("Preintervention"), drop = T)
names(NG.LF)
summary(NG.LF)

barplot(table(NG.LF$Occurrence)) # 260 absences and 932 presences

#Extracting the coord.ref for Nigeria 
sp <- SpatialPoints(NG.LF[,c('Longitude', 'Latitude')], proj4string = CRS("+proj=longlat +datum=WGS84"))
NG.LF <- SpatialPointsDataFrame(sp, NG.LF)
NG.LF

#This removes the preintervention column
NG.LF <- NG.LF[,1:3]

PCS <- Nigeria.ADM0@proj4string # extract Coord. System from Nigeria country map
NG.LF <- spTransform(NG.LF, PCS) # project the occurrence file
str(NG.LF)
class(NG.LF)
head(NG.LF[, 1:3])

#Plot the points using longitude and latitude
plot(Nigeria.ADM1)
points(NG.LF$Longitude, NG.LF$Latitude, pch=21, cex = 0.5, col="red") 

##Importing raster files 
# Here is how you return just the files names
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Covariates/Raster_covariates_initial") 

#Read the path to the files with a tif extension
raster.files <- list.files(pattern="*tif$", full.names=TRUE)
stopifnot(length(raster.files)>0)

raster.names <- list.files(pattern="*tif$", full.names=FALSE)
raster.names <- c(unlist(lapply(strsplit(raster.names,"[.]"), FUN=function(x) { x[1] })))   

# For loop assigning files names to individual raster objects
for (j in 1:length(raster.files)){
  assign(raster.names[j],raster(raster.files[j]))
}

rm(j)


# Rasterize Nigeria map to be use later with the imported covariates
raster.map <- raster(Nigeria.ADM0)
res(raster.map) <- res(AAcPrecip) # 1 km resolution
raster.map <- rasterize(Nigeria.ADM0, raster.map)
plot(raster.map)
plot(AAcPrecip)

# ~~~~~~~~~~~~ Loading pre-processed raster dataset from Covariate folder ~~~~~~~~~~~~~~~~~~~~~~~~ ##
#Prepare the stack raster object.
predictors <- stack(AAcPrecip, AMaxTemp, AMeanTemp, AMinTemp, AridityIndex, AnPET_AFRO, AvLST_AFRO, DQPrecip, WQPrecip,
                    MCQTemp, MWQTemp, EucDist_km_GLW3, EucDist_km_Rivers, FlowAccumulation, WI, EVI, Slope, 
                    SRTM1000_Africa, DMSP_2006_stable_lights, EucDist_2006_stable_lights_km, Sand_fraction,
                    Clay_fraction, Silt_fraction, pH)

names(predictors) <- c("Precipitation","Max Temp", "Mean Temp","Min Temp","AI","PET","LST","DQPrecip","WQPrecip",
                       "MCQTemp","MWQTemp","Distance Water Bodies","Distance Rivers","Flow Accumulation","Wettest Index",
                       "EVI","Slope","Elevation","Night-Light Emission","Distance Stable Lights","Sand","Clay","Silt",
                       "pH")

predictors <- mask(predictors, raster.map)

predictors <- stack(predictors) # Coerce to stack object

plot(predictors)

# Clean up the workspace from unnecessary raster objects (mostly raster layers loaded at the beginning)
rm(AAcPrecip, AMaxTemp, AMeanTemp, AMinTemp, AridityIndex, AnPET_AFRO, AvLST_AFRO, DQPrecip, WQPrecip,
   MCQTemp, MWQTemp, EucDist_km_GLW3, EucDist_km_Rivers, FlowAccumulation, WI, EVI, Slope, 
   SRTM1000_Africa, DMSP_2006_stable_lights, EucDist_2006_stable_lights_km, Sand_fraction,
   Clay_fraction, Silt_fraction, pH, ITN_2010_NG, ITN_2011_NG, ITN_2012_NG, ITN_2013_NG)


nlayers(predictors) # 28 potential predictors selected

gc()

#Explore layers included in the stack Object
for (i in 1:nlayers(predictors)){
  print(summary(predictors[[i]]))
}

rm(i)

## There are a few NAs cells in covariate layers that can trouble (i.e. dropping occurrence because void value) 
## as later during modelling. Run filter or focal statistics to fill NA cells. In case it's necessary
fill.na <- function(x, i=5) { # index [i] varies with the window set up for smoothing or correction (counting 1,0 down)
  if(is.na(x)[i] ) {
    return(round(mean(x, na.rm=TRUE),3))
  } else {
    return(round(x[i],3) )
  }
}  

predictors.corrected <- predictors
names(predictors.corrected)

for (j in 1:nlayers(predictors.corrected)) {
  predictors.corrected[[j]] <- focal(predictors.corrected[[j]], w = matrix(1,3,3), fun = fill.na, 
                                     pad = TRUE, na.rm = FALSE)
}

rm(j)

names(predictors.corrected) <- names(predictors)
names(predictors.corrected)
plot(predictors.corrected)

# Remove outliers in WQPrecip and DQPrecip raster layers
predictors.corrected[[8]] <- reclassify(predictors.corrected[[8]], cbind(65535, NA))
predictors.corrected[[9]] <- reclassify(predictors.corrected[[9]], cbind(65535, NA))

names(predictors.corrected) <- names(predictors)

predictors.corrected <- mask(predictors.corrected, raster.map)

predictors.corrected <- stack(predictors.corrected) # Coerce to stack object

predictors.corrected # 1,450,645 cells (1085 x 1337 matrix)

plot(predictors[[24]]) # displaying the mean temperature as background

plot(NG.LF, cex=0.5, col="red", add = TRUE) # adding in the LF occurrence points

#Checking collinearity in the covariates selected using functions on "usdm" package
#Based on Variance Inflation Factors (VIF). Calculate VIF for covariates
vif(predictors.corrected)

# identify variables using 'vifcor' or 'vifstep' that have little
# correlation based on correlation coefficients or VIFs- can be defined as
# an object and used to take an appropriate subset of variables
vifcor(predictors.corrected, th = 0.8)

# The above identifies variables with correlation coefficients < 0.80; can
# also be assigned to an object

# Alternatively, use the VIF to exclude correlated covariates.
# Create 'VIF' object listing variables that, together, have VIFs < 10 -
# this is a fairly standard threshold to use 

PredVars.NoCor <- vifstep(predictors.corrected, th = 10)

PredVars.NoCor

# 7 variables from the 24 input variables have collinearity problem: 
# Use the object created in the previous line to actually reduce the
# Predictor variables to be an uncorrelated set

Predictors.final <- exclude(predictors.corrected, PredVars.NoCor)

Predictors.final <- stack(Predictors.final)

names(Predictors.final) # 17 covariates left

#save rasters of environmental covariates 
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Covariates/Raster_covariates_final")

plot(Predictors.final$LST, main = "LST")
dev.print(file="LST.png", device=png, width = 900)

plot(Predictors.final$DQPrecip, main = "DQPrecip")
dev.print(file="DQPrecip.png", device=png, width = 900)

plot(Predictors.final$WQPrecip, main = "WQPrecip")
dev.print(file="WQPrecip.png", device=png, width = 900)

plot(Predictors.final$MCQTemp, main = "MCQTemp")
dev.print(file="MCQTemp.png", device=png, width = 900)

plot(Predictors.final$MWQTemp, main = "MWQTemp")
dev.print(file="MWQTemp.png", device=png, width = 900)

plot(Predictors.final$Distance.Water.Bodies, main = "Distance.Water.Bodies")
dev.print(file="Predictors.final.Distance.Water.Bodies.png", device=png, width = 900)

plot(Predictors.final$Distance.Rivers, main = "Distance.Rivers")
dev.print(file="Distance.Rivers.png", device=png, width = 900)

plot(Predictors.final$Flow.Accumulation, main = "Flow.Accumulation")
dev.print(file="Flow.Accumulation.png", device=png, width = 900)

plot(Predictors.final$Wettest.Index, main = "Wettest.Index")
dev.print(file="Wettest.Index.png", device=png, width = 900)

plot(Predictors.final$EVI, main = "EVI")
dev.print(file="EVI.png", device=png, width = 900)

plot(Predictors.final$Slope, main = "Slope")
dev.print(file="Slope.png", device=png, width = 900)

plot(Predictors.final$Elevation, main = "Elevation")
dev.print(file="Elevation.png", device=png, width = 900)

plot(Predictors.final$Night.Light.Emission, main = "Night.Light.Emission")
dev.print(file="Night.Light.Emission.png", device=png, width = 900)

plot(Predictors.final$Distance.Stable.Lights, main = "Distance.Stable.Lights")
dev.print(file="Distance.Stable.png", device=png, width = 900)

plot(Predictors.final$Clay, main = "Clay")
dev.print(file="Marginal_Plots_RF1.png", device=png, width = 900)

plot(Predictors.final$Silt, main = "Silt")
dev.print(file="Silt.png", device=png, width = 900)

plot(Predictors.final$pH, main = "pH")
dev.print(file="pH.png", device=png, width = 900)

# Clean up workspace

rm(fill.na, PredVars.NoCor)

#Extract covariates corresponding to coordinate location of LF occurence surveys
Covariates <- raster::extract(Predictors.final, NG.LF)
NG.LF.table <- cbind(as.data.frame(NG.LF),Covariates)
names(NG.LF.table)
head(NG.LF.table)
summary(NG.LF.table)

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/LF data")
write.csv(NG.LF.table, file = "LF.Env_Covariates.csv")

##Determine covariates to be included in the modelling process
## Using boosted regression trees to determine relative contribution of individual
##covariates to the occurence data. 
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/LF data")

install.packages("brt")
library(brt)
install.packages("dismo")
library(dismo)
install.packages("gbm")
library(gbm)


##Loading clean environmental variable and occurence data
LF.Env_Covariates <- read.csv("LF.Env_Covariates.csv")

##training data
LF.Env_Covariates.training <- gbm.step(data=LF.Env_Covariates,    
                                       gbm.x = c("LST", "DQPrecip", "WQPrecip", "MCQTemp", "MWQTemp", "Distance.Water.Bodies", "Distance.Rivers", 
                                                 "Flow.Accumulation", "Wettest.Index", "EVI", "Slope", "Elevation", "Night.Light.Emission", 
                                                 "Distance.Stable.Lights", "Clay", "Silt", "pH"),
                                       gbm.y = 3,
                                       family = "bernoulli",
                                       bag.fraction = 0.75,
                                       shrinkage = 0.005,
                                       n.trees = 40000,
                                       n.minobsinnode = 5,
                                       interaction.depth = 4,
                                       cv.folds = 10, 
                                       learning.rate = 0.001,
                                       tree.complexity = 10)

summary(LF.Env_Covariates.training)


##Simplifying the model
LF.Env_Covariates.training_simp <- gbm.simplify(LF.Env_Covariates.training, n.drops = 10)

##Now re-run the model indicating the number of dropped variables as in [[10]]

LF.Env_Covariates.training_dropped <- gbm.step(data=LF.Env_Covariates,    
                                               gbm.x = c("LST", "DQPrecip", "WQPrecip", "Wettest.Index", 
                                                         "Distance.Stable.Lights", "Slope", "Elevation"),
                                               gbm.y = 3,
                                               family = "bernoulli",
                                               bag.fraction = 0.75,
                                               shrinkage = 0.005,
                                               n.trees = 40000,
                                               n.minobsinnode = 5,
                                               interaction.depth = 4,
                                               cv.folds = 10,
                                               learning.rate = 0.001,
                                               tree.complexity = 10)

summary(LF.Env_Covariates.training_dropped)


##Plotting funtions and fitted values from the model
par(mfrow=c(1,1))
gbm.plot(LF.Env_Covariates.training_dropped, n.plots=7, write.title = F)

#Covariates that contributed <10% to the occurence data were excluded
names(Predictors.final)
Predictors.final2 <- dropLayer(Predictors.final, c(5, 4, 8, 17, 13, 15, 10, 8, 7, 16, 6))
names(Predictors.final2) #7 Covariates retained for fitting to ecological niche model. 


## ~~~~~~~~~~~~~~~~~~~~~~~~~Preparing validation dataset ~~~~~~~~~~~ #

#set working directory
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Spatially-stratified cross validation")

#load function which randomly selects validation dataset and grids up region of interest

do_boostrap <- function(dataset){
  #browser()
  idx <- unname(split(seq_len(nrow(dataset)), dataset$cell))
  pick <- sample(idx, size = length(idx), replace = TRUE)
  dataset[unlist(pick), ]
}

grid_up <- function(dataset, grid_size, rnd_dist){
  
  rd <- 0
  rd2 <- 0
  
  if (rnd_dist) {
    
    # draw random distance values 
    rd <- runif(n = 1, min = 0, max = grid_size)
    rd2 <- runif(n = 1, min = 0, max = grid_size)
    
  }
  
  # add rd to lat.grid and long.grid variables 
  dataset$lat.grid <- floor((dataset$Latitude - rd) / grid_size)    #replaced l for L in latitude
  dataset$long.grid <- floor((dataset$Longitude - rd2) / grid_size) #replaced l for L in longitude
  min.long <- min(dataset$long.grid)
  width.long <- max(dataset$long.grid) - min.long + 1
  min.lat <- min(dataset$lat.grid)
  
  dataset$cell <- (dataset$lat.grid - min.lat) * width.long + dataset$long.grid - min.long
  
  dataset
  
}



n_data <- 1192

grd_size <- 5 # decimal degrees (try 1 to 10 at least)

NG.LF.Spatial.Strat <- data.frame(id = seq_len(n_data), 
                                  Latitude = runif(n_data, 3, 14), 
                                  Longitude = runif(n_data, 2, 14), 
                                  Occurrence = runif(n_data, 0, 1))

NG.LF.gridded_dat <- grid_up(NG.LF.Spatial.Strat, grd_size, FALSE)

NG.LF.training_set <- do_boostrap(NG.LF.gridded_dat) # block bootstrap (one option)

NG.LF.train_indices <- unique(NG.LF.training_set$id)

NG.LF.valid_set <- NG.LF[-NG.LF.train_indices, ]

plot(Nigeria.ADM1)
points(NG.LF.valid_set, add = T, col = "red", cex = 0.5)

NG.LF.valid_set_occurrence <- NG.LF.valid_set[,3] #extracting only the column with Occurrence to feed into BIOMOD_FormatingData below


#Preparing predictors for validating dataset 

eval.covariates <- raster::extract(Predictors.final2, NG.LF.valid_set_occurrence)

eval.covariates <- as.data.frame(eval.covariates)

eval.covariates.table <- cbind(as.data.frame(NG.LF.valid_set_occurrence),eval.covariates)

eval.covariates.NG.LF <- SpatialPointsDataFrame(NG.LF.valid_set_occurrence, eval.covariates.table)

eval.covariates.NG.LF <- eval.covariates.NG.LF[,4:10]

## ~~~~~~~~~~~~~~~~~~~~~~~~~ Run Ecological Modelling using biomod2 package ~~~~~~~~~~~ #

## Compiling the different elements to be using when computing models
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/R_Scripts")

NG.LF.md <- NG.LF[,3]

NG.LF.Data.Formatted <- BIOMOD_FormatingData(resp.var = NG.LF.md,
                                             expl.var = Predictors.final2,
                                             resp.name = "Occurrence",
                                             eval.resp.var = NG.LF.valid_set_occurrence,
                                             eval.expl.var = eval.covariates.NG.LF)                                                                            


# Modelling. Building models based on the following algorithms: "GLM","GBM","ANN","SRE","MARS",
# "RF","MAXENT.Phillips"

MaxEnt <- "C:/Users/oae11/Desktop/NG_LF_ENM/R_scripts"

jar <- paste(MaxEnt, "maxent.jar", sep='/')

myBiomodOption <- BIOMOD_ModelingOptions(GBM = list(shrinkage = 0.005,
                                                    n.trees = 40000,
                                                    bag.fraction = 0.75,
                                                    n.minobsinnode = 5,
                                                    interaction.depth = 4,
                                                    cv.folds = 10,
                                                    learning.rate = 0.001, #Added learning rate 0.001
                                                    tree.complexity = 10), #Added tree complexity 10
                                         MAXENT.Phillips = list(MaxEnt.jar = jar))

# myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = jar)) if you want to use parameters by default, but always include point to MaxEnt app otherwise you won't get this SDM run.

system.time(NG.LF.model1 <- BIOMOD_Modeling(NG.LF.Data.Formatted,
                                            models = c("GLM","GBM","ANN","SRE","MARS",
                                                       "RF","MAXENT.Phillips"),
                                            models.options = myBiomodOption,
                                            NbRunEval = 100, # run at least 100 single models 
                                            DataSplit = 80, # spliting training/testing every round 
                                            Prevalence = NULL, 
                                            VarImport = 3,
                                            models.eval.meth = c("TSS","ROC","ACCURACY"),
                                            SaveObj = TRUE,
                                            rescal.all.models = FALSE,
                                            do.full.models = FALSE,
                                            modeling.id = paste("LF","ENM",sep =".")))
## system elapsed: 9431.04 (157.184 minutes)

NG.LF.model1 # none modells failing

##Extract evaluations from all models run
NG.LF.model1.Eval <- get_evaluations(NG.LF.model1)
dimnames(NG.LF.model1.Eval)

## Visualize performance of the models 
theme_set(theme_grey())
performance.models.1 <- models_scores_graph(NG.LF.model1, by = "models", metrics = c("ROC","TSS"),
                                            main = "Model performance comparison")

#making plot labellings bigger

performance.models.1 <- performance.models.1 + scale_color_manual(name="Model classes", 
                                                                  labels = c("ANN", 
                                                                             "GBM", 
                                                                             "GLM", 
                                                                             "MARS", 
                                                                             "MAXENT",
                                                                             "RF",
                                                                             "SRE"), 
                                                                  values = c("ANN"="blue", 
                                                                             "GBM"="red", 
                                                                             "GLM"="green", 
                                                                             "MARS"="brown", 
                                                                             "MAXENT.Phillips"="orange",
                                                                             "RF"="black",
                                                                             "SRE"="purple"))

performance.models.1 + theme_classic(base_size = 20)  + labs(x='Area under receiver operator curve', y='True skill statistics')

performance.models.1 + theme(legend.title = element_text(size = 20)) + labs(x='Area under receiver operator curve', y='True skill statistics', base_size = 36)

performance.models.1 + theme(legend.text = element_text(size = 18)) + labs(x='Area under receiver operator curve', y='True skill statistics')

performance.models.1 + theme(axis.title = element_text(size = 18)) + labs(x='Area under receiver operator curve', y='True skill statistics')

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Model_performace")

dev.print(file="Performance_models.png", device=png, width=1400)
dev.off()

performance.models.2 <- models_scores_graph(NG.LF.model1, by = "models", metrics = c("ROC","ACCURACY"))

dev.print(file="Performance_models2.png", device=png, width=900)
dev.off()

#For explanations of forecast scores - http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts
#Extract TSS/ROC/ACCURACY for test data (TSS: True Skill Statistic)
TSS.test.data <- as.data.frame(t(NG.LF.model1.Eval["TSS","Testing.data",,,]))
ROC.test.data <- as.data.frame(t(NG.LF.model1.Eval["ROC","Testing.data",,,]))
PCC.test.data <- as.data.frame(t(NG.LF.model1.Eval["ACCURACY","Testing.data",,,]))

summary(TSS.test.data)
summary(ROC.test.data)
summary(PCC.test.data)

#Extract TSS/ROC/ACCURACY for evaluation data (TSS: True Skill Statistic)
TSS.eval.data <- as.data.frame(t(NG.LF.model1.Eval["TSS","Evaluating.data",,,]))
ROC.eval.data <- as.data.frame(t(NG.LF.model1.Eval["ROC","Evaluating.data",,,]))
PCC.eval.data <- as.data.frame(t(NG.LF.model1.Eval["ACCURACY","Evaluating.data",,,]))

summary(TSS.eval.data)
summary(ROC.eval.data)
summary(PCC.eval.data)

TSS.test.Mean <- apply(TSS.test.data, 2, mean) 
TSS.test.Median <- apply(TSS.test.data, 2, median)
TSS.test.1Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
TSS.test.3Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.75), type = 7)) 

TSS.test <- t(data.frame(TSS.test.Mean, TSS.test.Median, TSS.test.1Q, TSS.test.3Q))

TSS.eval.Mean <- apply(TSS.eval.data, 2, mean) 
TSS.eval.Median <- apply(TSS.eval.data, 2, median)
TSS.eval.1Q <- apply(TSS.eval.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
TSS.eval.3Q <- apply(TSS.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

TSS.eval <- t(data.frame(TSS.eval.Mean, TSS.eval.Median, TSS.eval.1Q, TSS.eval.3Q))

ROC.test.Mean <- apply(ROC.test.data, 2, mean) 
ROC.test.Median <- apply(ROC.test.data, 2, median)
ROC.test.1Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
ROC.test.3Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.75), type = 7))

ROC.test <- t(data.frame(ROC.test.Mean, ROC.test.Median, ROC.test.1Q, ROC.test.3Q)) 

ROC.eval.Mean <- apply(ROC.eval.data, 2, mean) 
ROC.eval.Median <- apply(ROC.eval.data, 2, median)
ROC.eval.1Q <- apply(ROC.eval.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
ROC.eval.3Q <- apply(ROC.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

ROC.eval <- t(data.frame(ROC.eval.Mean, ROC.eval.Median, ROC.eval.1Q, ROC.eval.3Q)) 

PCC.test.Mean <- apply(PCC.test.data, 2, mean) 
PCC.test.Median <- apply(PCC.test.data, 2, median)
PCC.test.1Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.25), type = 7))
PCC.test.3Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.75), type = 7))

PCC.test <- t(data.frame(PCC.test.Mean, PCC.test.Median, PCC.test.1Q, PCC.test.3Q)) 

PCC.eval.Mean <- apply(PCC.eval.data, 2, mean) 
PCC.eval.Median <- apply(PCC.eval.data, 2, median)
PCC.eval.1Q <- apply(PCC.eval.data, 2, function (x) quantile(x, c(0.25), type = 7))
PCC.eval.3Q <- apply(PCC.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

PCC.eval <- t(data.frame(PCC.eval.Mean, PCC.eval.Median, PCC.eval.1Q, PCC.eval.3Q)) 

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Output tables")
write.csv(TSS.test, file = "TSS_test_NG_LF.csv")
write.csv(ROC.test, file = "ROC_test_NG_LF.csv")
write.csv(PCC.test, file = "PCC_test_NG_LF.csv")

write.csv(TSS.eval, file = "TSS_Eval_NG_LF.csv")
write.csv(ROC.eval, file = "ROC_Eval_NG_LF.csv")
write.csv(PCC.eval, file = "PCC_Eval_NG_LF.csv")

rm(TSS.test.data, ROC.test.data, PCC.test.data, TSS.test, TSS.eval.data, ROC.eval.data, PCC.eval.data, 
   TSS.eval, TSS.test.Mean, TSS.test.Median, TSS.test.1Q, TSS.test.3Q, TSS.eval.Mean, ROC.eval.Median,
   ROC.eval.1Q, ROC.eval.3Q, ROC.test, PCC.eval, PCC.eval.3Q, PCC.eval.1Q, PCC.eval.Median, PCC.eval.Mean, 
   PCC.test, PCC.test.3Q, PCC.test.1Q, PCC.test.Median, PCC.test.Mean,  ROC.eval, ROC.eval.Mean, ROC.test.1Q,
   ROC.test.3Q, ROC.test.Mean, ROC.test.Median, TSS.eval.1Q, TSS.eval.3Q, TSS.eval.Median)  

rm(TSS.Mean, TSS.Median, TSS.1Q, TSS.3Q, ROC.Mean, ROC.Median, ROC.1Q, ROC.3Q,
   PCC.Median, PCC.Mean, PCC.1Q, PCC.3Q)


setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Model_performace")

dev.print(file="Performance_models.png", device=png, width=1400)
dev.off()

performance.models.2 <- models_scores_graph(NG.LF.model1, by = "models", metrics = c("ROC","ACCURACY"))

dev.print(file="Performance_models2.png", device=png, width=900)
dev.off()

#For explanations of forecast scores - http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts
#Extract TSS/ROC/ACCURACY for test data (TSS: True Skill Statistic)
TSS.test.data <- as.data.frame(t(NG.LF.model1.Eval["TSS","Testing.data",,,]))
ROC.test.data <- as.data.frame(t(NG.LF.model1.Eval["ROC","Testing.data",,,]))
PCC.test.data <- as.data.frame(t(NG.LF.model1.Eval["ACCURACY","Testing.data",,,]))

summary(TSS.test.data)
summary(ROC.test.data)
summary(PCC.test.data)

#Extract TSS/ROC/ACCURACY for evaluation data (TSS: True Skill Statistic)
TSS.eval.data <- as.data.frame(t(NG.LF.model1.Eval["TSS","Evaluating.data",,,]))
ROC.eval.data <- as.data.frame(t(NG.LF.model1.Eval["ROC","Evaluating.data",,,]))
PCC.eval.data <- as.data.frame(t(NG.LF.model1.Eval["ACCURACY","Evaluating.data",,,]))

summary(TSS.eval.data)
summary(ROC.eval.data)
summary(PCC.eval.data)

TSS.test.Mean <- apply(TSS.test.data, 2, mean) 
TSS.test.Median <- apply(TSS.test.data, 2, median)
TSS.test.1Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
TSS.test.3Q <- apply(TSS.test.data, 2, function (x) quantile(x, c(0.75), type = 7)) 

TSS.test <- t(data.frame(TSS.test.Mean, TSS.test.Median, TSS.test.1Q, TSS.test.3Q))

TSS.eval.Mean <- apply(TSS.eval.data, 2, mean) 
TSS.eval.Median <- apply(TSS.eval.data, 2, median)
TSS.eval.1Q <- apply(TSS.eval.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
TSS.eval.3Q <- apply(TSS.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

TSS.eval <- t(data.frame(TSS.eval.Mean, TSS.eval.Median, TSS.eval.1Q, TSS.eval.3Q))

ROC.test.Mean <- apply(ROC.test.data, 2, mean) 
ROC.test.Median <- apply(ROC.test.data, 2, median)
ROC.test.1Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
ROC.test.3Q <- apply(ROC.test.data, 2, function (x) quantile(x, c(0.75), type = 7))

ROC.test <- t(data.frame(ROC.test.Mean, ROC.test.Median, ROC.test.1Q, ROC.test.3Q)) 

ROC.eval.Mean <- apply(ROC.eval.data, 2, mean) 
ROC.eval.Median <- apply(ROC.eval.data, 2, median)
ROC.eval.1Q <- apply(ROC.eval.data, 2, function (x) quantile(x, c(0.25), type = 7)) 
ROC.eval.3Q <- apply(ROC.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

ROC.eval <- t(data.frame(ROC.eval.Mean, ROC.eval.Median, ROC.eval.1Q, ROC.eval.3Q)) 

PCC.test.Mean <- apply(PCC.test.data, 2, mean) 
PCC.test.Median <- apply(PCC.test.data, 2, median)
PCC.test.1Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.25), type = 7))
PCC.test.3Q <- apply(PCC.test.data, 2, function (x) quantile(x, c(0.75), type = 7))

PCC.test <- t(data.frame(PCC.test.Mean, PCC.test.Median, PCC.test.1Q, PCC.test.3Q)) 

PCC.eval.Mean <- apply(PCC.eval.data, 2, mean) 
PCC.eval.Median <- apply(PCC.eval.data, 2, median)
PCC.eval.1Q <- apply(PCC.eval.data, 2, function (x) quantile(x, c(0.25), type = 7))
PCC.eval.3Q <- apply(PCC.eval.data, 2, function (x) quantile(x, c(0.75), type = 7))

PCC.eval <- t(data.frame(PCC.eval.Mean, PCC.eval.Median, PCC.eval.1Q, PCC.eval.3Q)) 

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Output tables")
write.csv(TSS.test, file = "TSS_test_NG_LF.csv")
write.csv(ROC.test, file = "ROC_test_NG_LF.csv")
write.csv(PCC.test, file = "PCC_test_NG_LF.csv")

write.csv(TSS.eval, file = "TSS_Eval_NG_LF.csv")
write.csv(ROC.eval, file = "ROC_Eval_NG_LF.csv")
write.csv(PCC.eval, file = "PCC_Eval_NG_LF.csv")

rm(TSS.test.data, ROC.test.data, PCC.test.data, TSS.test, TSS.eval.data, ROC.eval.data, PCC.eval.data, 
   TSS.eval, TSS.test.Mean, TSS.test.Median, TSS.test.1Q, TSS.test.3Q, TSS.eval.Mean, ROC.eval.Median,
   ROC.eval.1Q, ROC.eval.3Q, ROC.test, PCC.eval, PCC.eval.3Q, PCC.eval.1Q, PCC.eval.Median, PCC.eval.Mean, 
   PCC.test, PCC.test.3Q, PCC.test.1Q, PCC.test.Median, PCC.test.Mean,  ROC.eval, ROC.eval.Mean, ROC.test.1Q,
   ROC.test.3Q, ROC.test.Mean, ROC.test.Median, TSS.eval.1Q, TSS.eval.3Q, TSS.eval.Median)  

rm(TSS.Mean, TSS.Median, TSS.1Q, TSS.3Q, ROC.Mean, ROC.Median, ROC.1Q, ROC.3Q,
   PCC.Median, PCC.Mean, PCC.1Q, PCC.3Q)

##Get variable importance for the different models

NG.LF.model1.Var <- get_variables_importance(NG.LF.model1)
dimnames(NG.LF.model1.Var)

Variable.Contribution <- t(apply(NG.LF.model1.Var, 2, rowMeans))

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Output tables")
write.csv(Variable.Contribution, file = "Variable_Contribution_LF.csv")


setwd("C:/Users/oae11/Desktop/NG_LF_ENM/R_Scripts")

##Explore variable response curves for GBM and RF models
NF.LF.GBM <- BIOMOD_LoadModels(NG.LF.model1, models = "GBM") 

GBMsPlot2D <- response.plot2(models = NF.LF.GBM, Data = get_formal_data(NG.LF.model1,"expl.var"), 
                             show.variables = get_formal_data(NG.LF.model1, "expl.var.names"), 
                             do.bivariate = FALSE, fixed.var.metric = "median", col = c("blue", "red"),legend = TRUE, 
                             data_species = get_formal_data(NG.LF.model1, "resp.var"))

NF.LF.RF <- BIOMOD_LoadModels(NG.LF.model1, models = "RF")

RFsPlot2D <- response.plot2(models = NF.LF.RF, Data = get_formal_data(NG.LF.model1,"expl.var"), 
                            show.variables = get_formal_data(NG.LF.model1, "expl.var.names"), 
                            do.bivariate = FALSE, fixed.var.metric = "median", col = c("blue", "red"),legend = TRUE, 
                            data_species = get_formal_data(NG.LF.model1, "resp.var"))

rm(list = ls(pattern = "Occurrence_AllData")) 

Selected.models <- c(NF.LF.GBM, NF.LF.RF)

#Ensemble Modelling. Come up with a consensus model and ensemble for every type of model
NG.LF.model1.Ensemble <- BIOMOD_EnsembleModeling(modeling.output = NG.LF.model1,
                                                 chosen.models = Selected.models,
                                                 em.by = "all",
                                                 eval.metric = c("ROC"),
                                                 eval.metric.quality.threshold = c(0.85),  
                                                 prob.mean = T,
                                                 prob.cv = T,
                                                 prob.ci = T,
                                                 prob.ci.alpha = 0.05,
                                                 prob.median = T,
                                                 committee.averaging = T,
                                                 prob.mean.weight = T,
                                                 prob.mean.weight.decay = "proportional")

# chosen.models must match the selection for projection. Models can be selected by passing on a
# vector with models' name for those that may want to be considered for modelling.
# For instance, BRT (GBM) and RF are the only models which exceed ROC value of 0.725, so ensure to 
#include only ensemble which fulfill this criteria based on these two modelling approach.

## eval.metric.quality.threshold: evaluation threshold to exclude individual models for being assembled.
## 0.725 means models ROC < 0.725 will be disregarded, and not included in the ensemble.

NG.LF.model1.Ensemble

#Extract evaluation of the ensemble model

Ensemble.Algorithms <- c("Mean","Coef.Var","InfCI","SupCI","Median","MCAvg","WeightedMean")

NG.LF.Ensemble.Eval <- get_evaluations(NG.LF.model1.Ensemble)

names(NG.LF.Ensemble.Eval)

NG.LF.Ensemble.TSS <- t(sapply(NG.LF.Ensemble.Eval, "[",2,1:5))
rownames(NG.LF.Ensemble.TSS) <- Ensemble.Algorithms

NG.LF.Ensemble.KAPPA <- t(sapply(NG.LF.Ensemble.Eval, "[",1,1:5))
rownames(NG.LF.Ensemble.KAPPA) <- Ensemble.Algorithms

NG.LF.Ensemble.ROC <- t(sapply(NG.LF.Ensemble.Eval, "[",3,1:5))
rownames(NG.LF.Ensemble.ROC) <- Ensemble.Algorithms

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Output tables")
write.csv(NG.LF.Ensemble.TSS, file = "NG_LF_Ensemble_TSS.csv")
write.csv(NG.LF.Ensemble.KAPPA, file = "NG_LF_Ensemble_KAPPA.csv")
write.csv(NG.LF.Ensemble.ROC, file = "NG_LF_Ensemble_ROC.csv")

#setwd(path_original)
setwd("C:/Users/oae11/Desktop/NG_LF_ENM/R_Scripts")

#Final projection (extrapolation) based on the best fitted model. 
## Space projection for the algorithms selected
gc()

NG.LF.Proj <- BIOMOD_Projection(modeling.output = NG.LF.model1,
                                new.env = Predictors.final2,
                                proj.name = "Final",
                                selected.models = Selected.models,
                                binary.meth = "ROC",
                                compress = "gzip",
                                build.clamping.mask = FALSE,
                                output.format = ".grd")

## There is no point of projecting for models which have not been ensembled. Hence,
# project for "Selected models" (i.e. BRT+RF) and disregard the rest.

NG.LF.Proj
list.files("Occurrence/proj_Final")

NG.LF.projections <- get_predictions(NG.LF.Proj) # Set up a raster stack with the predictions

#Ensemble forescasting combining projections based on models ensemble rules defined at the ensemble modelling step
NG.LF.EF <- BIOMOD_EnsembleForecasting(EM.output = NG.LF.model1.Ensemble,
                                       projection.output = NG.LF.Proj)

NG.LF.Ensemble <- get_predictions(NG.LF.EF)

names(NG.LF.Ensemble)

# Extract best performing models: mean of probabilities with CI and models committee averaging
# probability raster are provided scale 0 to 1000, so rescale from 0 to 1
# This is made because float and double raster are heavier (and take up more disk space)

NG.LF.mean <- raster(NG.LF.Ensemble, layer= 1)/1000 
NG.LF.LB <- raster(NG.LF.Ensemble, layer= 3)/1000
NG.LF.UB <- raster(NG.LF.Ensemble, layer = 4)/1000
NG.LF.median <- raster(NG.LF.Ensemble, layer = 5)/1000
NG.LF.MCA <- raster(NG.LF.Ensemble, layer = 6)/1000
NG.LF.WM <- raster(NG.LF.Ensemble, layer = 7)/1000

par(mfrow=c(1,1))
plot(NG.LF.LB)
plot(Nigeria.ADM1, add= TRUE)
plot(NG.LF.median)
plot(Nigeria.ADM1, add= TRUE)
plot(NG.LF.UB)
plot(Nigeria.ADM1, add= TRUE)


writeRaster(NG.LF.mean,filename = "NG_LF_Mean.tif", format= "GTiff", overwrite=TRUE)
writeRaster(NG.LF.LB,filename = "NG_LF_LB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(NG.LF.UB,filename = "NG_LF_UB.tif", format= "GTiff", overwrite=TRUE)
writeRaster(NG.LF.median,filename = "NG_LF_median.tif", format= "GTiff", overwrite=TRUE)
writeRaster(NG.LF.MCA,filename = "NG_LF_MCA.tif", format= "GTiff", overwrite=TRUE)
writeRaster(NG.LF.WM,filename = "NG_LF_WM.tif", format= "GTiff", overwrite=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##Creating marginal effect plots for covariates modelled using BRT
theme_set(theme_grey())


LST <- gather(GBMsPlot2D$LST, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM) 

LST.plot.GBM <- ggplot(LST, aes(x= LST, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Land Surface Temperature (?C) RC - 11.29%") + ylab("Response Curve") +
  theme_bw()

LST.plot.GBM <- LST.plot.GBM + theme_classic(base_size = 12)

LST.plot.GBM <- LST.plot.GBM + labs(x='Degrees celsius (RC 11.29%)', title='Land Surface Temperature')


DQPrecip <- gather(GBMsPlot2D$DQPrecip, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM)

DQPrecip.plot.GBM <- ggplot(DQPrecip, aes(x= DQPrecip, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Driest quarter precipitation (mm) RC - 28.99%") +   ylab("Response Curve") +
  theme_bw()

DQPrecip.plot.GBM <- DQPrecip.plot.GBM + theme_classic(base_size = 12)

DQPrecip.plot.GBM <- DQPrecip.plot.GBM + labs(x='Millimetres (RC 28.99%)', title='Driest quarter precipitation')


WQPrecip <- gather(GBMsPlot2D$WQPrecip, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM) 

WQPrecip.plot.GBM <- ggplot(WQPrecip, aes(x= WQPrecip, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Wettest quarter precipitation (mm) RC - 18.13%") +  ylab("Response Curve") +
  theme_bw()

WQPrecip.plot.GBM <- WQPrecip.plot.GBM + theme_classic(base_size = 12)

WQPrecip.plot.GBM <- WQPrecip.plot.GBM + labs(x='Millimetres (RC 18.13%)', title='Wettest quarter precipitation')


Distance.Stable.Lights <- gather(GBMsPlot2D$Distance.Stable.Lights, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM)

Distance.Stable.Lights.plot.GBM <- ggplot(Distance.Stable.Lights, aes(x= Distance.Stable.Lights, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Distance to stable lights (km) RC - 6.61%") +  ylab("Response Curve") +
  theme_bw()

Distance.Stable.Lights.plot.GBM <- Distance.Stable.Lights.plot.GBM + theme_classic(base_size = 12)

Distance.Stable.Lights.plot.GBM <- Distance.Stable.Lights.plot.GBM + labs(x='Kilometres (RC 6.61%)', title='Distance to stable lights')


Wettest.Index <- gather(GBMsPlot2D$Wettest.Index, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM)

Wettest.Index.plot.GBM <- ggplot(Wettest.Index, aes(x= Wettest.Index, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Wetness index RC - 13.21%") +   ylab("Response Curve") +
  theme_bw()

Wettest.Index.plot.GBM <- Wettest.Index.plot.GBM + theme_classic(base_size = 12)

Wettest.Index.plot.GBM <- Wettest.Index.plot.GBM + labs(x='Wetness index (RC 13.21%)', title='Wetness index')


Slope <- gather(GBMsPlot2D$Slope, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM)

Slope.plot.GBM <- ggplot(Slope, aes(x= Slope, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Terrain slope (Degrees) RC - 7.99%") +   ylab("Response Curve") +
  theme_bw()

Slope.plot.GBM <- Slope.plot.GBM + theme_classic(base_size = 12)

Slope.plot.GBM <- Slope.plot.GBM + labs(x='Degree (RC 7.99%)', title='Terrain slope')


Elevation <- gather(GBMsPlot2D$Elevation, Model, Marginal_Effect, Occurrence_AllData_RUN1_GBM:Occurrence_AllData_RUN100_GBM)

Elevation.plot.GBM <- ggplot(Elevation, aes(x= Elevation, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Elevation (m) RC - 13.79%") +   ylab("Response Curve") +
  theme_bw()

Elevation.plot.GBM <- Elevation.plot.GBM + theme_classic(base_size = 12)

Elevation.plot.GBM <- Elevation.plot.GBM + labs(x='Metres above sea level (RC 13.79%)', title='Elevation')


# Produce and export combined plots

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Marginal_effects_plots")

plot_grid(DQPrecip.plot.GBM, WQPrecip.plot.GBM, Elevation.plot.GBM, Wettest.Index.plot.GBM,  
          LST.plot.GBM, Slope.plot.GBM, Distance.Stable.Lights.plot.GBM, labels = c("a","b","c","d","e","f","g"),
          ncol = 3, nrow = 3, label_size = 12)


dev.print(file="Marginal_Plots_GBM.png", device=png, width = 900)
dev.off()



##Creating marginal effect plots for covariates modelled using RF
theme_set(theme_grey())

LST <- gather(RFsPlot2D$LST, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF) #RUN cganged to 10 from 100

LST.plot.RF <- ggplot(LST, aes(x= LST, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Land Surface Temperature (?C) RC- 15.05%") + ylab("Response Curve") +
  theme_bw()

LST.plot.RF <- LST.plot.RF + theme_classic(base_size = 12)

LST.plot.RF <- LST.plot.RF + labs(x='Degrees celsius (RC 15.05%)', title='Land Surface Temperature')


DQPrecip <- gather(RFsPlot2D$DQPrecip, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF)

DQPrecip.plot.RF <- ggplot(DQPrecip, aes(x= DQPrecip, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Driest quarter precipitation (mm) RC - 27.95%") +   ylab("Response Curve") +
  theme_bw()

DQPrecip.plot.RF <- DQPrecip.plot.RF + theme_classic(base_size = 12)

DQPrecip.plot.RF <- DQPrecip.plot.RF + labs(x='Millimetres (RC 27.95%)', title='Driest quarter precipitation')


WQPrecip <- gather(RFsPlot2D$WQPrecip, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF) #RUN cganged to 10 from 100

WQPrecip.plot.RF <- ggplot(WQPrecip, aes(x= WQPrecip, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Wettest quarter precipitation (mm) RC - 16.47%") +  ylab("Response Curve") +
  theme_bw()

WQPrecip.plot.RF <- WQPrecip.plot.RF + theme_classic(base_size = 12)

WQPrecip.plot.RF <- WQPrecip.plot.RF + labs(x='Millimetres (RC 16.47%)', title='Wettest quarter precipitation')


Distance.Stable.Lights <- gather(RFsPlot2D$Distance.Stable.Lights, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF)

Distance.Stable.Lights.plot.RF <- ggplot(Distance.Stable.Lights, aes(x= Distance.Stable.Lights, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Distance to stable lights (km) RC - 8.94%") +  ylab("Response Curve") +
  theme_bw()

Distance.Stable.Lights.plot.RF <- Distance.Stable.Lights.plot.RF + theme_classic(base_size = 12)

Distance.Stable.Lights.plot.RF <- Distance.Stable.Lights.plot.RF + labs(x='Kilometres (RC 8.94%)', title='Distance to stable lights')


Wettest.Index <- gather(RFsPlot2D$Wettest.Index, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF)

Wettest.Index.plot.RF <- ggplot(Wettest.Index, aes(x= Wettest.Index, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Wetness Index RC - 9.52%") +   ylab("Response Curve") +
  theme_bw()

Wettest.Index.plot.RF <- Wettest.Index.plot.RF + theme_classic(base_size = 12)

Wettest.Index.plot.RF <- Wettest.Index.plot.RF + labs(x='Wetness Index (RC 9.52%)', title='Wetness Index')


Slope <- gather(RFsPlot2D$Slope, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF)

Slope.plot.RF <- ggplot(Slope, aes(x= Slope, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Terrain slope (Degrees) RC - 7.10%") +   ylab("Response Curve") +
  theme_bw()

Slope.plot.RF <- Slope.plot.RF + theme_classic(base_size = 12)

Slope.plot.RF <- Slope.plot.RF + labs(x='Degree (RC 7.10%)', title='Terrain slope')


Elevation <- gather(RFsPlot2D$Elevation, Model, Marginal_Effect, Occurrence_AllData_RUN1_RF:Occurrence_AllData_RUN100_RF)

Elevation.plot.RF <- ggplot(Elevation, aes(x= Elevation, y= Marginal_Effect)) +
  scale_y_continuous(limits = c(0.5,1)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal,
               fun.args = list(conf.int= 0.95), fill= "grey") +
  stat_summary(geom = "line", fun.y = mean, colour = "blue", size = 1) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size= 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.position = "none") +
  xlab("Elevation (m) RC - 14.96%") +   ylab("Response Curve") +
  theme_bw()

Elevation.plot.RF <- Elevation.plot.RF + theme_classic(base_size = 12)

Elevation.plot.RF <- Elevation.plot.RF + labs(x='Metres above sea level (RC 14.96%)', title='Elevation')


# Produce and export combined plots

setwd("C:/Users/oae11/Desktop/NG_LF_ENM/Outputs and covariates/Marginal_effects_plots")

plot_grid(DQPrecip.plot.RF, WQPrecip.plot.RF, LST.plot.RF, Elevation.plot.RF, Wettest.Index.plot.RF,   
          Distance.Stable.Lights.plot.RF, Slope.plot.RF, labels = c("a","b","c","d","e","f","g"),
          ncol = 3, nrow = 3, label_size = 12)


dev.print(file="Marginal_Plots_RF.png", device=png, width = 900)
dev.off()

