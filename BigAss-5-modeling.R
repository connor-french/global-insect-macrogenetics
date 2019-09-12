##### MARS Modeling
library(raster)
library(tidymodels)
library(tidyverse)
library(earth)
library(pdp)

###Specify the results folder
#if you want to work out of a folder from an earlier date, replace this string with the date
todays_date <- Sys.Date()

#folder for the entire project's output to go into
todays_results <- paste0("results_", todays_date)

###Start modeling
#read in genetic data
sum_df <- read_csv("sum_df_2019-06-14.csv")

###read in environmental rasters. I have to hard code this because I can't host them on the github and they need to be downloaded
raster_folder <-"/Users/connorfrench/Dropbox/Old_Mac/climate-data/chelsa_10min/10min" 

#list all of the tif rasters
f <- list.files(raster_folder, 
                pattern = "tif$",
                full.names = TRUE)

#read in the files as a stack
envs <- stack(f)

#resample the env't to 100 km cells. This takes a while for global datasets
sa_clim_1d <- envs %>% 
  resample_equal_area(km = 100)

envs <- stack("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/BigAss-bird-phylogeography/ndvi_dhi_combined-v6/ndvi_dhi_combo_downscaled_1-degree.tif") %>% 
  stack(envs)



#extract environmental data for modeling
coordinates(sum_df) <- ~longitude+latitude
env_data <- raster::extract(envs, sum_df) %>% 
  as_tibble()

sum_df_combo <- sum_df %>% 
  as.data.frame(xy = TRUE) %>% 
  as_tibble() %>% 
  bind_cols(env_data) %>% 
  filter(median.pi.avg < 0.05, mean.pi.avg < 0.05, longitude < 0) %>% 
  rename(dhi_1 = ndvi_dhi_combo_downscaled_1.degree.1,
         dhi_2 = ndvi_dhi_combo_downscaled_1.degree.2,
         dhi_3 = ndvi_dhi_combo_downscaled_1.degree.3) %>% 
  drop_na()
  

#set up recipe
model_recipe <- recipe(hill.one.avg ~ 
                         latitude + 
                         longitude + 
                         bio1 + 
                         bio2 + 
                         bio3 + 
                         bio4 + 
                         bio5 + 
                         bio6 + 
                         bio7 +
                         bio8 +
                         bio9 +
                         bio10 +
                         bio11 +
                         bio12 +
                         bio13 +
                         bio14 +
                         bio15 +
                         bio16 +
                         bio17 +
                         bio18 +
                         bio19 +
                         dhi_1 +
                         dhi_2 +
                         dhi_3, 
                       data = sum_df_combo
                       ) 
  

#prepare the recipe so it can be applied to other data sets (usually do this with training data that is going to be applied to other testing sets) 
#since we're just doing the training and validation steps (no testing) and aren't transforming the variables at all, this doesn't really matter, but I'm doing it for good practice
prepped_recipe <- prep(model_recipe, sum_df_combo)

#apply the recipe to a data set to create a design matrix
sum_df_processed <- bake(prepped_recipe, sum_df_combo)

#create a feature plot to visualize relationships between variables
caret::featurePlot(select(sum_df_processed, latitude, longitude, starts_with("bio"), contains("dhi")), sum_df_processed$hill.one.avg)


#run mars algorithm
mars_model <- earth(
  x = select(sum_df_processed, latitude, longitude, starts_with("bio"), contains("dhi")),
  y = sum_df_processed$hill.one.avg,
  degree = 1
)

plot(mars_model)



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
)
  
model_fit_plot <- model_fit_df %>% 
  mutate(dataset = fct_reorder(dataset, rsq)) %>% 
  ggplot(aes(x = dataset, y = rsq)) +
  labs(x = "", y = "Rsq (GRsq)", title = "Model fit for geographic subsets") +
  geom_col(fill = "#440154FF") +
  geom_col(aes(x = dataset, y = grsq), fill = "#20A486FF") +
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



