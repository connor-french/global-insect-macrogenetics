MARS modeling of genetic diversity ~ environment. First model is a MARS model using generalized cross-validation
Load packages
```{r include = FALSE}
library(tidymodels)
library(earth)
library(caret)
library(ggrepel)
```


REDO THE ANALYSES WITH A RECIPE
```{r}
mars.data <- pi_env %>% 
  raster::as.data.frame(xy = TRUE) %>% 
  filter(!duplicated(cells)) %>% 
  dplyr::select(hill.one.avg, contains("bio_10m")) %>% 
  na.omit()


#split into training and testing data
mars.split <- initial_split(mars.data, prop = 4/5)

###########create training data set##############
mars.training <- training(mars.split)

#select bioclims as features
mars.training.features <- mars.training %>% dplyr::select(contains("bio_10m")) 

#select mean pi as response
mars.training.response <- mars.training %>% dplyr::select(hill.one.avg) %>% as.matrix() %>% as.numeric()

#list of feature dataframe and response matrix
mars.training <- list(features = mars.training.features, response =  mars.training.response) 

############create testing data set############
mars.test <- testing(mars.split) 

#select bioclims as features
mars.test.features <- mars.test %>% dplyr::select(contains("bio_10m")) 

#select mean pi as response
mars.test.response <- mars.test %>% dplyr::select(hill.one.avg) %>% as.matrix() %>% as.numeric()

#list of feature dataframe and response matrix
mars.test <- list(features = mars.test.features, response =  mars.test.response) 



featurePlot(mars.training$features, mars.training$response)

#gcv model
mars.fit.gcv <- earth(mars.training$features, mars.training$response)
mars.pred.gcv <- predict(mars.fit.gcv, mars.test$features)

mars.plot.data <- tibble(pred = as.numeric(mars.pred.gcv), obs = as.numeric(mars.test$response))

ggplot(mars.plot.data, aes(x = obs, y = pred)) +
  geom_point() +
  geom_abline(alpha = 0.5)

mars.fit.gcv

```


```{r}
#5-fold cross validation
mars.grid <- expand.grid(degree = 1:2, nprune = seq(1, 19, by = 1))

#LOOCV
ctrl <- trainControl(
  method = "LOOCV", 
  number = 10,
  # Save the assessment predictions from the best model
  savePredictions = "final",
  # Log the progress of the tuning process
  verboseIter = FALSE
)

# Using the same seed to obtain the same 
# resamples.
set.seed(223489)
mars_mod <- train(
  mars.training$features, 
  mars.training$response,
  method = "earth",
  tuneGrid = mars.grid,
  trControl = ctrl
)

ggplot(mars_mod) + theme(legend.position = "top")
```


```{r}
mars_imp <- varImp(mars_mod)
ggplot(mars_imp, top = 20) + xlab("")

ggplot(mars_mod$pred, aes(x = pred, y = obs)) +
  geom_point() +
  geom_abline()
```


prediction performance
```{r}
mars.test.cv <- mars.test$features %>% 
  as_tibble() %>% 
  mutate(pred = predict(mars_mod, .$features))# %>% 
summarise(rmse = rmse_vec(., truth = obs, estimate = pred), rsq = rsq_vec(., truth = obs, estimate = pred))

```
