список packages

library(sp)
library(rgdal)
library(raster)
library(dynatopmodel)
library(spatialEco)
library(gstat)
library(automap)
library(tidyverse)
library(readr)
library(sf)

points <- read.csv('xrf_points_with_nickel.csv')
dem <- readGDAL('predictors/dsm_p_m_c_1.5.tif')

subset <- subset(points, Set == "Test")
validation <- subset(points, Set == "Validation")

subset$Ni_exp_sqrt <- sqrt(subset$Ni_express)

coordinates(subset) <- ~x_utm + y_utm

proj4string(subset) <-  CRS("+init=epsg:32636")
plot(subset)

predictors <- raster(dem)
predictors$dem <- raster(dem)
predictors$slope <- terrain(raster(dem), opt = 'slope', unit = 'degrees')
predictors$aspect <- terrain(raster(dem), opt = 'aspect', unit = 'degrees')
predictors$flowdir <- terrain(raster(dem), opt = 'flowdir')
predictors$flowacc <- upslope.area(raster(dem), log=F, atb=F, deg=0.1, fill.sinks=F)
predictors$twi <- log(predictors$flowacc / tan(predictors$slope / 180))
predictors$curvature <- curvature(raster(dem), type = "total")

predictors$soils <- raster(readGDAL('predictors/mon_soils_raster_1.5m.tif'))

preds <- dropLayer(predictors, 1)
predictors <- preds
proj4string(predictors) <-  CRS("+init=epsg:32636")

plot(predictors)

pred <- as(predictors, 'SpatialGridDataFrame')
pred$soils <- as.factor(pred$soils)

sample <- raster::extract(predictors, subset, df=TRUE)
sample$soils <- as.factor(sample$soils)

write.csv2(subset, "test_subset.csv")
write.csv2(sample, "test_sample.csv")

subset$slope <- sample$slope
subset$aspect <- sample$aspect
subset$dem <- sample$dem
subset$curvature <- sample$curvature
subset$twi <- sample$twi
subset$flowdir <- sample$flowdir
subset$flowacc <- sample$flowacc
subset$soils <- sample$soils

write.csv2(subset, "test_subset2.csv")
write.csv2(sample, "test_sample2.csv")

a3 = data.frame(subset)
a4 = na.exclude(a3) 
subset <- subset(a4)
coordinates(subset) <- ~x_utm + y_utm
proj4string(subset) <-  CRS("+init=epsg:32636")

head(subset)

l_fit <- lm(Ni_exp_sqrt ~ slope + aspect + dem + curvature + twi + flowdir + flowacc + soils, subset)

optimal_fit <- step(l_fit, direction = 'backward')

summary(optimal_fit)

vgmlm <- variogram(optimal_fit$residuals ~ 1, subset)


vgm_lm <- vgm(nugget = 0.4, psill = 0.35, range = 60, model = "Exp")

plot(vgmlm, vgm_lm)

rk <- krige(optimal_fit$call$formula, subset, pred, vgm_lm)

krig <- autoKrige(Ni_exp_sqrt~1, subset, pred)

plot(krig)
plot(krig$krige_output$var1.pred)

pred$PredNiLmRK <- (rk$var1.pred)^2 

pred$OrKrig <- (krig$krige_output$var1.pred)^2

predicted <- stack(pred)
plot(predicted)

coordinates(validation) <- ~x_utm + y_utm
proj4string(validation) <-  CRS("+init=epsg:32636")

sampleVAL <- raster::extract(predicted, validation, df = T)

validation$Ni_predicted_lm_rk <- sampleVAL$PredNiLmRK
validation$Ni_predicted_or_kg <- sampleVAL$OrKrig

plot(validation$Ni_express, validation$Ni_predicted_lm_rk)
plot(validation$Ni_express, validation$Ni_predicted_or_kg)

cor(validation$Ni_express, validation$Ni_predicted_lm_rk)
cor(validation$Ni_express, validation$Ni_predicted_or_kg)
