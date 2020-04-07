##### Predictive modeling interactive script
##### This is an interactive script. 
library(raster)
library(tidyverse)
library(earth)
library(caret)
library(pdp)
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


#### repeat the process for all resolutions of climate data (100 km, 200 km, 300 km)
full_df <- extract_env(env_raster_dir = "data/climate/rasters_100km",
                       gen_data_path = "data/genetics/sum_df_100.csv")

# change the filename for each resolution
write_csv(full_df, "data/genetics/full_df_100.csv")
 
# read in full df if you've already extracted the environment values
full_df <- read_csv("data/genetics/full_df_100.csv")

#### create feature plots to visualize relationships between variables
#### have to do one for each data set or the plot gets way too busy

# chelsa current
caret::featurePlot(dplyr::select(full_df, 
                          latitude, 
                          longitude, 
                          starts_with("chelsa")), 
                   full_df$hill.three.avg)

# dhi
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("dhi")), 
                   full_df$hill.one.avg)

# heterogeneity
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("hetero")), 
                   full_df$hill.one.avg)

# late Holocene
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_2")), 
                   full_df$hill.one.avg)

# mid Holocene
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_3")), 
                   full_df$hill.one.avg)

# younger dryas stadial
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_4")), 
                   full_df$hill.one.avg)

# bolling-allerod
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_5")), 
                   full_df$hill.one.avg)

# heinrich stadial
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_6")), 
                   full_df$hill.one.avg)

# last glacial maximum
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("pc_7")), 
                   full_df$hill.one.avg)

# soil
caret::featurePlot(dplyr::select(full_df, 
                                 starts_with("soil")), 
                   full_df$hill.one.avg)


responses <- full_df %>%
  na.omit() %>% 
  dplyr::select(cell,
         latitude,
         longitude,
         contains("avg"),
         contains("sd"),
         contains("skew"))

predictors <- full_df %>% 
  na.omit() %>% 
  dplyr::select(
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  )

cor_df <- full_df %>% 
  na.omit() %>% 
  correlate(method = "spearman") %>% 
  stretch()

response_half <- cor_df %>% 
  na.omit() %>% 
  filter(
    r > 0.3,
    str_detect(x, "cell") |
    str_detect(x, "latitude") |
    str_detect(x, "longitude") |
    str_detect(x, "avg") |
    str_detect(x, "sd") |
    str_detect(x, "skew")
    )



# none of the genetic response variables have a correlation > 0.3 with a non-genetic variable
response_half %>% View()

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

#####################################

### Hill 1 model (use all vars) #####

#####################################

# first hill number analysis
full_df_hill <- full_df %>%
  dplyr::select(
    hill.one.avg,
    latitude,
    longitude,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  ) %>%
  na.omit()

#### MARS model

# set up tuning grid. We're exploring additive and interaction terms ("degree"), 
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 50, by = 2))


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

mars_model <- mod

# average RMSE plot
ggplot(mars_model) + theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred <- mars_model$pred

ggplot(mars_model_pred, aes(x = pred, y = obs)) +
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



####### Random forest

# set up tuning grid. We're exploring additive and interaction terms ("degree"), 
rf_grid <- expand.grid(.mtry = 1:sqrt(ncol(full_df_hill)))

set.seed(8305445)
rf_model <- caret::train(
  x = select(full_df_hill, -hill.one.avg),
  y = full_df_hill$hill.one.avg,
  method = "rf",
  tuneGrid = rf_grid,
  trControl = ctrl,
  ntree = 1000
)

rf_model$results



#####################################

### Avg. pi model (use all vars) #####

#####################################

# first hill number analysis
full_df_pi <- full_df %>%
  select(
    mean.pi.avg,
    latitude,
    longitude,
    starts_with("chelsa"),
    starts_with("dhi"),
    starts_with("hetero"),
    starts_with("pc"),
    starts_with("soil")
  ) %>%
  na.omit()

#### MARS model

# set up tuning grid. We're exploring additive and interaction terms ("degree"), 
mars_grid <- expand.grid(degree = 1:2, nprune = seq(2, 50, by = 2))

set.seed(3453465)
mars_model_pi <- caret::train(
  x = select(full_df_pi, -mean.pi.avg),
  y = responses,
  method = "earth",
  tuneGrid = mars_grid,
  trControl = ctrl
)

# average RMSE plot
ggplot(mars_model_pi) + theme(legend.position = "top")

# obs vs predicted plot
mars_model_pred_pi <- mars_model_pi$pred

ggplot(mars_model_pred_pi, aes(x = pred, y = obs)) +
  geom_abline(col = "green", alpha = .5) +
  geom_point(alpha = .8) +
  geom_smooth(
    se = FALSE,
    col = "red",
    lty = 2,
    lwd = .25,
    alpha = .5
  ) +
  xlim(0, 10)

####### Random forest

# set up tuning grid. We're exploring additive and interaction terms ("degree"), 
set.seed(833472)
rf_model_pi <- caret::train(
  x = select(full_df_pi, -mean.pi.avg),
  y = full_df_pi$mean.pi.avg,
  method = "rf",
  tuneGrid = rf_grid,
  trControl = ctrl,
  ntree = 1000
)

rf_model_pi$results

######################################################

################### Extra Plots ######################

######################################################
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








