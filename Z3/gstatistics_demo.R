library(sp) #
library(rgdal) # чтение форматов геоданных
library(raster) # работа с растровыми данными
library(dynatopmodel) # морфометрический анализ
library(spatialEco) # обработка пространственных данных
library(gstat) # геостатистика
library(automap) # кригинг

points <- read.csv('xrf_points_demo.csv')
dem <- readGDAL('predictors/dsm_p_m_c_1.5.tif')

# Разбиваем выборку на 2 части (~75% - тестовые данные, ~25% - данные для валидации)

subset <- subset(points, Set == "Test")
validation <- subset(points, Set == "Validation")

# Данные для линейных регрессионных моделей должны быть распределены нормально
# В противном случае нужно данные трансформировать

subset$Cu_exp_sqrt <- sqrt(subset$Cu_express)

# Конвертируем таблицу в набор пространственных данных

coordinates(subset) <- ~x_utm + y_utm

# Очень важно присвоить правильную проекцию
# https://www.earthdatascience.org/courses/earth-analytics/spatial-data-r/understand-epsg-wkt-and-other-crs-definition-file-types/

proj4string(subset) <-  CRS("+init=epsg:32636")
plot(subset)

# Готовим растровые данные - независимые переменные для моделирования распределения концентрации меди в почве

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

# Для регрессионного кригинга нам нужно конвертировать наш набор растров в формат 'SpatialGridDataFrame'

pred <- as(predictors, 'SpatialGridDataFrame')
pred$soils <- as.factor(pred$soils)


# Для каждой точки мы теперь можем "вытащить" значения каждого параметра в соотвествующей ячейки растра

sample <- raster::extract(predictors, subset, df=TRUE)
sample$soils <- as.factor(sample$soils)
subset$slope <- sample$slope
subset$aspect <- sample$aspect
subset$dem <- sample$dem
subset$curvature <- sample$curvature
subset$twi <- sample$twi
subset$flowdir <- sample$flowdir
subset$flowacc <- sample$flowacc
subset$soils <- sample$soils

head(subset)

# Построим нашу первую линейную регрессионную модель

l_fit <- lm(Cu_exp_sqrt ~ slope + aspect + dem + curvature + twi + flowdir + flowacc + soils, subset)

optimal_fit <- step(l_fit, direction = 'backward')

summary(optimal_fit)

# Создаем вариограмму остатков

vgmlm <- variogram(optimal_fit$residuals ~ 1, subset)

# Также необходимо создать оптимальную модель вариограммы

vgm_lm <- vgm(nugget = 0.07, psill = 0.085, range = 60, model = "Sph")

plot(vgmlm, vgm_lm)

# Регрессионный кригинг

rk <- krige(optimal_fit$call$formula, subset, pred, vgm_lm)
# это пока модель

# Заодно сделаем обычный кригинг
krig <- autoKrige(Cu_exp_sqrt~1, subset, pred)

plot(krig)
plot(krig$krige_output$var1.pred)

# Теперь построим карты концентрации меди в почве на основе полученных моделей

pred$PredCuLmRK <- (rk$var1.pred)^2 # почему мы возводим тут результат в квадрат?

pred$OrKrig <- (krig$krige_output$var1.pred)^2

# Конвертируем обратно в набор растров

predicted <- stack(pred)
plot(predicted)

# Подготовим данные для валидации

coordinates(validation) <- ~x_utm + y_utm
proj4string(validation) <-  CRS("+init=epsg:32636")

sampleVAL <- raster::extract(predicted, validation, df = T)

validation$Cu_predicted_lm_rk <- sampleVAL$PredCuLmRK
validation$Cu_predicted_or_kg <- sampleVAL$OrKrig

# Проанализируем корреляцию предсказанных значений концентрации меди и данных полевых измерений

plot(validation$Cu_express, validation$Cu_predicted_lm_rk)
plot(validation$Cu_express, validation$Cu_predicted_or_kg)

cor(validation$Cu_express, validation$Cu_predicted_lm_rk)
cor(validation$Cu_express, validation$Cu_predicted_or_kg)

# Можно также посчитать среднюю ошибку между двуми столбцами
# Какой метод лучше моделирует зависимую переменную, почему ?

