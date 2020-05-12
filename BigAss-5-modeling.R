##### Predictive modeling interactive script
##### This is an interactive script. 
library(raster)
library(earth)
library(caret)
library(pdp)
library(ClustOfVar)
library(here)
library(tidyverse)
library(corrr)
### Helper functions

# extract environment for lat-longs and add to genetics df
extract_env <- function(env_raster_dir, gen_data_path) {
  gen_data <- readr::read_csv(gen_data_path)
  sp::coordinates(gen_data) <- ~longitude+latitude
  
  env_stack <- list.files(env_raster_dir, pattern = ".tif$", full.names = TRUE) %>% 
    raster::stack()
  
  env_data <- raster::extract(env_stack, gen_data) %>% 
    tibble::as_tibble()
  
  all_data <- gen_data %>% 
    raster::as.data.frame(xy = TRUE) %>% 
    tibble::as_tibble() %>% 
    dplyr::bind_cols(env_data)
  
  return(all_data)
}


### extract environmental values for each resolution ------
# resolution you want to extract at
resolution <- 100


full_df_path <- here("data", 
                     "genetics", 
                     paste0("full_df_", resolution, ".csv"))


# if the extracted data frame is already there, then just read in the data frame
if (file.exists(full_df_path)) {
  full_df <- read_csv(full_df_path)
} else {
  env_raster_dir <- here("data", "climate", paste0("rasters_", resolution, "km"))
  gen_data_path <- here("data", "genetics", paste0("sum_df_", resolution, ".csv"))
  full_df <- extract_env(env_raster_dir = env_raster_dir,
                         gen_data_path = gen_data_path)
  
  # write to file
  write_csv(full_df, full_df_path)
}

# we want the absolute value of latitude to be included in the model (distance from equator) 
full_df <- full_df %>% mutate(
  abs_lat = abs(latitude),
  abs_long = abs(longitude)
)

# check out NA values. Most come from the earthenv heterogeneity layers. All in all, not too many NAs (23 total)
full_df %>% map_dbl(function(x)
  sum(length(which(is.na(x))))) %>% 
  enframe(name = "variable", value = "na_count") %>% 
  filter(na_count > 0)

### Determine data subsets to keep (remove colinearity and maximize information) -----
## Modifying the Crowther lab's script from their Nematode paper. (Find script here: https://gitlab.ethz.ch/devinrouth/crowther_lab_nematodes/blob/master/Nematode_Geospatial/Nematode_ClustOfVar.ipynb)

# select only predictor variables
predictors <- full_df %>% 
  na.omit() %>% 
  dplyr::select(
    abs_lat,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )


# Run the clustering function and plot the variables
predictor_tree <- hclustvar(X.quanti=as.matrix(predictors))
plot(predictor_tree)


# Enter desired number of clusters to generate (in a vector)
num_clusters <-  c(5,10,25, 50)

# Prep an empy list to store the variables from the clustering
vars_list <-  list()


# Loop the clustering across the desired number of clusters
for (nC in num_clusters) {
  # Set the random number seed
  set.seed(347)
  
  # Examine the top variables in each cluster
  part <- cutreevar(predictor_tree, nC)
  
  # Instantiate an empty list in which to store the variables
  mid_variable_list <-  c()
  
  # Extract only the names of the variables (i.e., unlist the values in their parent data structure)
  names_list_with_nulls <- lapply(part$var, rownames)
  names_list_without_nulls <-
    Filter(Negate(function(x)
      is.null(unlist(x))), names_list_with_nulls)
  
  # Ascertain the top variables from the KMeans Clustering
  # and store it in a list
  for (j in c(1:length(names_list_without_nulls))) {
    mid_variable_list[j] <-  names_list_without_nulls[[j]][1]
  }
  
  # Prep the variable list for input into the model calls
  vars_to_include <- list((mid_variable_list))
  vars_list[match(nC, num_clusters)] <- vars_to_include
  
}

predictors_5 <- predictors %>% select(vars_list[[1]])

predictors_10 <- predictors %>% select(vars_list[[2]])

predictors_25 <- predictors %>% select(vars_list[[3]])
 
predictors_50 <- predictors %>% select(vars_list[[4]])


responses <- full_df %>%
  na.omit() %>% 
  dplyr::select(cell,
                abs_lat,
                abs_long,
                contains("avg"),
                contains("sd"),
                contains("skew"))


######### Set up models ##########

# we're performing 10-fold cross validation for all models
ctrl <- trainControl(
  method = "cv", 
  number = 10,
  # Save the assessment predictions from the best model
  savePredictions = "final",
  # Log the progress of the tuning process
  verboseIter = TRUE,
  # want to use parallelization 
  allowParallel = TRUE 
)

# set up tuning grid. We're exploring additive and interaction terms ("degree"), 
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 50, by = 2))


### Model performance across responses -----

mars_model_fun <- function(x) {
  caret::train(
    x = predictors,
    y = responses[[x]],
    method = "earth",
    tuneGrid = mars_grid,
    trControl = ctrl
  )
  
}

for (name in colnames(responses)) {
  mod <- mars_model_fun(name)
  save(mod, file = paste0("model_output_11-12/", name, ".rda"))
}

load(here("archive", "model_output_11-12", "hill.eight.avg.rda"))

mod_pred <- mod$pred

ggplot(mod_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### Hill 1 model  -----

### Full Model -----
# need to change variable name to mars_model_full
full_mars_model <- caret::train(
  x = predictors,
  y = responses$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(full_mars_model) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_full <- full_mars_model$pred

ggplot(mars_model_pred_full, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )


### Five predictor model -----
mars_model_5 <- caret::train(
  x = predictors_5,
  y = responses$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_5) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_5 <- mars_model_5$pred

ggplot(mars_model_pred_5, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )


### 10 predictor model -----
mars_model_10 <- caret::train(
  x = predictors_10,
  y = responses$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_10) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_10 <- mars_model_10$pred

ggplot(mars_model_pred_10, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### 25 predictor model -----
mars_model_25 <- caret::train(
  x = predictors_25,
  y = responses$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_25) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_25 <- mars_model_25$pred

ggplot(mars_model_pred_25, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )


### 50 predictor model -----
mars_model_50 <- caret::train(
  x = predictors_50,
  y = responses$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_50) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_50 <- mars_model_50$pred

ggplot(mars_model_pred_50, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )



### Random forest -----

### Full Model -----
# set up tuning grid. 
rf_grid <- expand.grid(.mtry = 1:sqrt(ncol(predictors)))

set.seed(8445)
rf_model <- caret::train(
  x = predictors,
  y = responses$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid,
  trControl = ctrl,
  ntree = 1000
)

rf_model$results

# average RMSE plot
ggplot(rf_model) + 
  theme(legend.position = "top")

# obs vs predicted plot
rf_model_pred <- rf_model$pred

ggplot(rf_model_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### 5 predictor model -----

rf_grid_5 <- expand.grid(.mtry = 1:ncol(predictors_5))

set.seed(348)
rf_model_5 <- caret::train(
  x = predictors_5,
  y = responses$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid_5,
  trControl = ctrl,
  ntree = 1000
)

rf_model_5$results

# average RMSE plot
ggplot(rf_model_5) + 
  theme(legend.position = "top")

# obs vs predicted plot
rf_model_5_pred <- rf_model_5$pred

ggplot(rf_model_5_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### 10 predictor model -----

rf_grid_10 <- expand.grid(.mtry = 1:ncol(predictors_10))

set.seed(38)
rf_model_10 <- caret::train(
  x = predictors_10,
  y = responses$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid_10,
  trControl = ctrl,
  ntree = 1000
)

rf_model_10$results

# average RMSE plot
ggplot(rf_model_10) + 
  theme(legend.position = "top")

# obs vs predicted plot
rf_model_10_pred <- rf_model_10$pred

ggplot(rf_model_10_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### 25 predictor model -----

rf_grid_25 <- expand.grid(.mtry = 1:ncol(predictors_25))

set.seed(25)
rf_model_25 <- caret::train(
  x = predictors_25,
  y = responses$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid_25,
  trControl = ctrl,
  ntree = 1000
)

rf_model_25$results

# average RMSE plot
ggplot(rf_model_25) + 
  theme(legend.position = "top")

# obs vs predicted plot
rf_model_25_pred <- rf_model_25$pred

ggplot(rf_model_25_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### 50 predictor model -----

rf_grid_50 <- expand.grid(.mtry = 1:15)

set.seed(50)
rf_model_50 <- caret::train(
  x = predictors_50,
  y = responses$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid_50,
  trControl = ctrl,
  ntree = 1000
)

rf_model_50$results

# average RMSE plot
ggplot(rf_model_50) + 
  theme(legend.position = "top")

# obs vs predicted plot
rf_model_50_pred <- rf_model_50$pred

ggplot(rf_model_50_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### Regional models ----

### Northern hemisphere -----
# N = 488
predictors_north <- full_df %>% 
  na.omit() %>% 
  filter(latitude > 0) %>% 
  dplyr::select(
    abs_lat,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )

responses_north <- full_df %>%
  filter(latitude > 0) %>% 
  na.omit() %>% 
  dplyr::select(cell,
                abs_lat,
                abs_long,
                contains("avg"),
                contains("sd"),
                contains("skew"))

mars_model_north <- caret::train(
  x = predictors_north,
  y = responses_north$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_north) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_north_pred <- mars_model_north$pred

ggplot(mars_model_north_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### Southern hemisphere -----
# N = 151
predictors_south <- full_df %>% 
  na.omit() %>% 
  filter(latitude < 0) %>% 
  dplyr::select(
    abs_lat,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )

responses_south <- full_df %>%
  filter(latitude < 0) %>% 
  na.omit() %>% 
  dplyr::select(cell,
                abs_lat,
                abs_long,
                contains("avg"),
                contains("sd"),
                contains("skew"))

mars_model_south <- caret::train(
  x = predictors_south,
  y = responses_south$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_south) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_south_pred <- mars_model_south$pred

ggplot(mars_model_south_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )


### Western hemisphere -----
# N = 328
predictors_west <- full_df %>% 
  na.omit() %>% 
  filter(longitude < 0) %>% 
  dplyr::select(
    abs_lat,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )

responses_west <- full_df %>%
  filter(longitude < 0) %>% 
  na.omit() %>% 
  dplyr::select(cell,
                abs_lat,
                abs_long,
                contains("avg"),
                contains("sd"),
                contains("skew"))

mars_model_west <- caret::train(
  x = predictors_west,
  y = responses_west$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_west) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_west_pred <- mars_model_west$pred

ggplot(mars_model_west_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )

### Eastern hemisphere -----
# N = 311
predictors_east <- full_df %>% 
  na.omit() %>% 
  filter(longitude > 0) %>% 
  dplyr::select(
    abs_lat,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )

responses_east <- full_df %>%
  filter(longitude > 0) %>% 
  na.omit() %>% 
  dplyr::select(cell,
                abs_lat,
                abs_long,
                contains("avg"),
                contains("sd"),
                contains("skew"))

mars_model_east <- caret::train(
  x = predictors_east,
  y = responses_east$hill.one.avg,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_east) + 
  theme(legend.position = "top")

# obs vs predicted plot
mars_model_east_pred <- mars_model_east$pred

ggplot(mars_model_east_pred, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_hex(alpha = .8) +
  scale_fill_viridis_c() +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  )









### Extra Plots -----

var_importance <- tibble(
  variable = dimnames(evimp(mars_model))[[1]],
  num_subsets = c(22, 21, 18, 15, 14, 12, 10, 9, 6, 5)
)

#variable importance barplot
var_imp_plot <- var_importance %>% 
  mutate(variable = fct_reorder(variable, num_subsets)) %>% 
  ggplot(aes(x = variable, y = num_subsets)) +
  labs(x = "Environmental variable", y = "Number of Subsets", title = "Variable importance for global dataset") +
  geom_col(fill = "#440154FF") +
  coord_flip() +
  theme(text = element_text(size=35, color = "#333333"),
        panel.background = element_rect(fill = "transparent", color = "#333333"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        legend.text=element_text(color = "#333333"),
        legend.title = element_text(color = "#333333"),
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg)
        axis.text.x = element_text(color = "#333333"),
        axis.text.y = element_text(color = "#333333")
  )

dhi_plot <- ggplot(sum_df_processed, aes(x = dhi_1, y = hill.one.avg)) +
  geom_point(color = "gold") +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  labs(x = "DHI 3: Seasonality", y = "Avg Hill One") +
  theme(text = element_text(size=35, color = "white"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        legend.text=element_text(color = "white"),
        legend.title = element_text(color = "white"),
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", color = "white"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg)
        axis.text.x = element_text(color = "white"),
        axis.text.y = element_text(color = "white")
  )


#model fit across hemispheres
model_fit_df <- tibble(
  dataset = c("Global", "North_Hemisphere", "South_Hemisphere", "East_Hemisphere", "West_Hemisphere"),
  rsq = c(0.2422463, 0.1824755, 0.5103241, 0.2856198, 0.2800988),
  grsq = c(0.1831, 0.1366418, 0.3169007, 0.2056086, 0.2157808),
  n = c(1321, 893, 301, 620, 574)
) %>% 
  pivot_longer(cols = c(rsq, grsq), names_to = "metric", values_to = "value")

model_fit_plot <- model_fit_df %>% 
  mutate(dataset = fct_reorder(dataset, value)) %>% 
  ggplot(aes(x = dataset, y = value, fill = metric)) +
  labs(x = "", y = "Rsq (GRsq)", title = "Model fit for geographic subsets") +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#440154FF", "#20A486FF")) +
  coord_flip() +
  theme(text = element_text(size=20, color = "white"),
        panel.background = element_rect(fill = "transparent", color = "white"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        legend.text=element_text(color = "white"),
        legend.title = element_text(color = "white"),
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg)
        axis.text.x = element_text(color = "white"),
        axis.text.y = element_text(color = "white")
  )


ggsave(var_imp_plot, 
        filename = "/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/Conferences/Evolution-2019/poster/var_imp_plot.png", 
        bg= "transparent",
        device = "png")

ggsave(dhi_plot, 
       filename = "/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/Conferences/Evolution-2019/poster/dhi_plot.png", 
       bg= "transparent",
       device = "png")

ggsave(model_fit_plot, 
       filename = "/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/Conferences/Evolution-2019/poster/model_fit_plot.png", 
       bg= "transparent",
       device = "png")








