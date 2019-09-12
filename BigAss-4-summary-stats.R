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

### combine into a single raster stack
clim <- stack(chelsa_agg, dhi_agg, ghh_agg)


# filter sequences based on: 
# pi (0.05, to make it similar to threshold used for BOLD, with some leeway),
# number of sequences >= 5 (based on literature recommendations)
# percent missing < 0.25 (general rule of thumb)

pi_df_one <- 
  pi_df_one %>% 
  filter(avg_pi <= 0.05, num_seqs >= 5, perc_missing < 0.25) %>%
  group_by(cell) %>% 
  filter(n() >= 10)

#plot density plots of avg pi and sd pi
avg_pi_density <- 
  pi_df_one %>% 
  ggplot(aes(x = avg_pi, y = ..density..)) +
  geom_histogram(fill = "darkgreen", color = "black", bins = 50) +
  geom_density(fill = "gray", alpha = 0.2) +
  labs(x = "Average pi", y = "Density", title = "Average pi per species per cell") +
  theme_minimal() 

sd_pi_density <- 
  pi_df_one %>% 
  ggplot(aes(x = sd_pi, y = ..density..)) +
  geom_histogram(fill = "darkgreen", color = "black", bins = 50) +
  geom_density(fill = "gray", alpha = 0.2) +
  labs(x = "SD pi", y = "Density", title = "SD pi per species per cell") +
  theme_minimal() 

#make plot list
d_plot_grob <- list(avg_pi_density, sd_pi_density) %>% lapply(ggplotGrob)

#save plots to output
ggsave(paste0(plots, "/pi_density_plots.pdf"), gridExtra::marrangeGrob(grobs = d_plot_grob, nrow = 2, ncol = 1))

#replicate the sampling 1000 times. 
sum_list <- purrr::rerun(1000, pi_summary_fun(pi_df_one, n_species = 10)) 

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


#' 
#' 

