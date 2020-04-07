###### Plot sequences
#library(furrr)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)
library(moments)
library(here)
library(tidyverse)

#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
source("R/pad_short_seqs_function.R") #pad short sequences with Ns
source("R/fasta_write_function.R") #write fasta files 
source("R/calc_gen_dist.R") #align seqs and calculate mean pi per species per cell
source("R/hill_num_calc.R") #calculate hill number of pi. Note* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.
source("R/pi_summary_function.R")
source("R/sumstat_plot_function.R")

read_env <- function(env_data_path) {
  env_data <- list.files(env_data_path, pattern = ".tif$", full.names = TRUE) %>% 
    raster::stack()
  
  return(env_data)
}


env_data <- read_env("data/climate/rasters_100km")

pi_df_raw <- read_csv("data/genetics/pi_df_one_100.csv")

test_path <- "../results_2019-07-07_100/output/test_nuc.csv"
test_nuc <- 
  data.table::fread(test_path)
data.table::setnames(test_nuc, "cells", "cell")
# filter sequences based on: 
# number of sequences >= 5 (based on literature recommendations)
# percent missing < 0.5 (general rule of thumb)
# number of BINs >= 10

pi_df <- 
  pi_df_raw %>% 
  filter(num_seqs >= 5, perc_missing < 0.5) %>% 
  group_by(cell) %>% 
  filter(n() >= 10)


#replicate the sampling 100 times. 
sum_list <- purrr::rerun(100, pi_summary_fun(pi_df, n_species = 10)) 

#summarise these samples into a final df
sum_df <- bind_rows(sum_list) %>%
  group_by(cell) %>% 
  dplyr::summarise(
    median.pi.avg = mean(median.pi), 
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
    hill.three.avg = mean(hill.three),
    hill.three.sd = sd(hill.three),
    hill.four.avg = mean(hill.four),
    hill.four.sd = sd(hill.four),
    hill.five.avg = mean(hill.five),
    hill.five.sd = sd(hill.five),
    hill.six.avg = mean(hill.six),
    hill.six.sd = sd(hill.six),
    hill.seven.avg = mean(hill.seven),
    hill.seven.sd = sd(hill.seven),
    hill.eight.avg = mean(hill.eight),
    hill.eight.sd = mean(hill.eight),
    shannon.avg = mean(shannon), 
    shannon.sd = sd(shannon)
  )

sum_df <- read_csv("data/genetics/sum_df_100.csv")
#add lat longs back to df
sum_df_latlongs <- test_nuc %>% 
  as_tibble() %>% 
  dplyr::filter(!duplicated(cell)) %>% 
  dplyr::select(cell, latitude = latitude, longitude = longitude) %>% 
  right_join(sum_df, by = "cell") 


write_csv(sum_df_latlongs, "sum_df_100.csv")
#sum_df <- fread("pi-summary-insects-10.csv")




####### PLOTS ################

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
  pi_raster_one <- rasterize(sum_df_latlongs, clim, fun = "first", field = stat)
  writeRaster(pi_raster_one, filename = paste0("rasters_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".tif"), format = "GTiff")
  #leaflet plot. 
  l <- leaflet() %>% 
    addTiles() %>%
    addRasterImage(pi_raster_one, colors = pal_one,  opacity = 0.8) %>%
    addLegend(pal = pal_one, values = values(pi_raster_one))
  #write to plots folder
  mapview::mapshot(l, selfcontained = FALSE, url = paste0(getwd(), "/rasters_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".html"))
}


# plot density plots of avg pi and sd pi
avg_pi_density <- 
  pi_df %>% 
  ggplot(aes(x = avg_pi, y = ..density..)) +
  geom_histogram(fill = "darkgreen", color = "black", bins = 50) +
  geom_density(fill = "gray", alpha = 0.2) +
  labs(x = "Average pi", y = "Density", title = "Average pi per species per cell") +
  theme_minimal() 

#make plot list
d_plot_grob <- list(avg_pi_density, sd_pi_density) %>% lapply(ggplotGrob)

#save plots to output
ggsave(paste0(plots, "/pi_density_plots.pdf"), gridExtra::marrangeGrob(grobs = d_plot_grob, nrow = 2, ncol = 1))


