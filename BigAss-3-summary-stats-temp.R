###### Plot sequence
library(data.table)
library(tidyverse)
library(tidylog)
#library(furrr)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)
library(moments)

#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
source("R/pad_short_seqs_function.R") #pad short sequences with Ns
source("R/fasta_write_function.R") #write fasta files 
source("R/calc_gen_dist.R") #align seqs and calculate mean pi per species per cell
source("R/hill_num_calc.R") #calculate hill number of pi. Note* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.
source("R/pi_summary_function.R")
source("R/sumstat_plot_function.R")

###Specify the results folder
#if you want to work out of a folder from an earlier date, replace this string with the date
todays_date <- Sys.Date()

#folder for the entire project's output to go into
todays_results <- paste0("results_", todays_date)

###Begin summary stats
#read in sequence statistics dataframe. replace path if you're conducting this analysis on a different day than the first two steps
pi_df_file <- paste0("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/BigAss-bird-phylogeography/BigAss-phylogeography/pi_df_one_2019-05-26.csv")
pi_df_one <- fread(pi_df_file) %>% 
  mutate(cells = cell)

pi_df_one %>% 
  #filter(avg_pi > 0.01) %>% 
  group_by(cells) %>% 
  count(sort = TRUE)

#read in test_nuc dataframe.
test_nuc_path <- paste0("test_nuc_2019-05-21.csv")
test_nuc <- fread(test_nuc_path)

#read in environmental data for plotting
envs <- raster::getData(name = "worldclim", var = "bio", res = 10)

#read in bioclim rasters, downscale the resolution to 1 degree, and crop them
sa_clim_1d <- envs %>% #stack(f) %>% 
  #  crop(bounds) %>% 
  aggregate(fact = 6)


#plot density plots of avg pi and sd pi
avg_pi_density <- pi_df_one %>% 
  ggplot(aes(x = avg_pi)) +
  geom_density(fill = "darkgreen", alpha = 0.7) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7, bins = 50) +
  labs(x = "Average pi", y = "Frequency", title = "Average pi per species per cell") +
  theme_minimal() 

sd_pi_density <- pi_df_one %>% 
  ggplot(aes(x = sd_pi)) +
  geom_density(fill = "darkgreen", alpha = 0.7) +
  geom_histogram(fill = "darkgreen", color = "black") +
  labs(x = "SD pi", y = "Frequency", title = "SD pi per species per cell") +
  theme_minimal() 

#make plot list
d_plot_grob <- list(avg_pi_density, sd_pi_density) %>% lapply(ggplotGrob)

#save plots to output
ggsave(paste0("plots_", Sys.Date(), "/pi_density_plots_", Sys.Date(), ".pdf"), gridExtra::marrangeGrob(grobs = d_plot_grob, nrow = 2, ncol = 1))

#replicate the sampling 1000 times. 
sum_list <- replicate(1000, pi_summary_fun(pi_df_one, 10), simplify = FALSE) 

#summarise these samples into a final df
sum_df <- bind_rows(sum_list) %>%
  group_by(cells) %>% 
  dplyr::summarise(median.pi.avg = mean(median.pi), 
                   median.pi.sd = sd(median.pi), 
                   median.pi.skew = skewness(median.pi),
                   mean.pi.avg = mean(mean.pi),  
                   mean.pi.sd = sd(mean.pi), 
                   mean.pi.skew = skewness(mean.pi),
                   sd.pi.avg = mean(sd.pi), 
                   sd.pi.sd = sd(sd.pi), 
                   hill.zero.avg = mean(hill.zero), 
                   hill.zero.sd = sd(hill.zero), 
                   hill.one.avg = mean(hill.one), 
                   hill.one.sd = sd(hill.one), 
                   hill.two.avg = mean(hill.two), 
                   hill.two.sd = sd(hill.two), 
                   shannon.avg = mean(shannon), 
                   shannon.sd = sd(shannon))

#add lat longs back to df
sum_df_latlongs <- test_nuc %>% 
  as_tibble() %>% 
  filter(!duplicated(cells)) %>% 
  dplyr::select(cells, latitude = Lat, longitude = Long) %>% 
  right_join(sum_df, by = "cells") %>% 
  filter(mean.pi.avg < 0.05) #%>% #remove species-level divergence
  

sum_df_latlongs %>% 
  #filter(mean.pi.skew < .015) %>% 
  ggplot(aes(x = abs(latitude), y = hill.one.avg)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_minimal()


fwrite(sum_df_latlongs, file = paste0("sum_df_", Sys.Date(), ".csv"))
#sum_df <- fread("pi-summary-insects-10.csv")


#write df to file
fwrite(pi_env, paste0("pi_env_", Sys.Date(), ".csv"))



#extract the relevant column names (any name that contains a summary of pi, the hill number, or shannon entropy)
gen_sum_vec <- grep("pi|hill|shannon", colnames(sum_df_latlongs), value = TRUE)
gen_sum_vec <- set_names(gen_sum_vec)

#loop through columns and plot histograms
pi_sum_stat_plots <- map(gen_sum_vec, ~sumstat_plot_fun(data = sum_df_latlongs, .x)) %>% 
  map(ggplotGrob)

#save plots to output
ggsave(paste0("pi_sumstat_plots_", Sys.Date(), ".pdf"), gridExtra::marrangeGrob(grobs = pi_sum_stat_plots, nrow = 2, ncol = 1))


#convert back to a spatial data frame
coordinates(sum_df_latlongs) <- ~longitude+latitude

####plot maps of the genetic summary stats on the landscape
#color palette for leaflet plot
pal_one <- colorNumeric(palette = "viridis", domain = NULL, na.color = "#00000000")


#plot genetic summary stats on a map and output the raster
dir.create(paste0("rasters_", Sys.Date()))
for(stat in gen_sum_vec){
  #create raster of genetic summary stat per cell
  pi_raster_one <- rasterize(sum_df_latlongs, sa_clim_1d, fun = "first", field = stat)
  writeRaster(pi_raster_one, filename = paste0("rasters_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".tif"), format = "GTiff")
  #leaflet plot. 
  l <- leaflet() %>% 
    addTiles() %>%
    addRasterImage(pi_raster_one, colors = pal_one,  opacity = 0.8) %>%
    addLegend(pal = pal_one, values = values(pi_raster_one))
  #write to plots folder
  mapview::mapshot(l, selfcontained = FALSE, url = paste0(getwd(), "/rasters_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".html"))
}
 

#' 
#' 


summary(lm(sum_df_processed$hill.one.avg~sum_df_processed$bio2))
