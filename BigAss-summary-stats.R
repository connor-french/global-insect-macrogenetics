######

#plot density plots of avg pi and sd pi
avg_pi_density <- pi_df_one %>% 
  ggplot(aes(x = avg_pi)) +
  geom_density(fill = "darkgreen", alpha = 0.7) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7) +
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
sum_list<- replicate(1000, pi_summary_fun(pi_df_one, 10), simplify = FALSE) 

#summarise these samples into a final df
sum_df <- bind_rows(sum_list) %>%
  group_by(cells) %>% 
  dplyr::summarise(median.pi.avg = mean(median.pi), median.pi.sd = sd(median.pi), mean.pi.avg = mean(mean.pi),  mean.pi.sd = sd(mean.pi), sd.pi.avg = mean(sd.pi), sd.pi.sd = sd(sd.pi), hill.zero.avg = mean(hill.zero), hill.zero.sd = sd(hill.zero), hill.one.avg = mean(hill.one), hill.one.sd = sd(hill.one), hill.two.avg = mean(hill.two), hill.two.sd = sd(hill.two), shannon.avg = mean(shannon), shannon.sd = sd(shannon))


fwrite(sum_df, file = paste0("output_", Sys.Date(), "/sum_df_", Sys.Date(), ".csv"))
#sum_df <- fread("pi-summary-insects-10.csv")

#make new data frame including the pi values.
pi_env <- merge(test_nuc, sum_df, by = "cells") %>% 
  raster::as.data.frame(xy = TRUE) %>% #have to convert to data frame to omit NAs
  dplyr::select(recordID, bin_uri, cells, markercode, Long, Lat, contains("pi"), contains("hill"), contains("shannon"), contains("bio")) %>%
  #filter(mean.pi.avg < 0.035) %>% #cells largest pi values contain "species" which actually contain different species
  #na.omit() %>%
  drop.levels()

#write df to file
fwrite(pi_env, paste0("output_", Sys.Date(), "/pi_env_", Sys.Date(), ".csv"))



#extract the relevant column names (any name that contains a summary of pi, the hill number, or shannon entropy)
gen_sum_vec <- grep("pi|hill|shannon", colnames(pi_env), value = TRUE)
gen_sum_vec <- set_names(gen_sum_vec)

#loop through columns and plot histograms
pi_sum_stat_plots <- map(gen_sum_vec, ~sumstat_plot_fun(data = pi_env, .x)) %>% 
  map(ggplotGrob)

#save plots to output
ggsave(paste0("plots_", Sys.Date(), "/pi_sumstat_plots_", Sys.Date(), ".pdf"), gridExtra::marrangeGrob(grobs = pi_sum_stat_plots, nrow = 2, ncol = 1))


#convert back to a spatial data frame
coordinates(pi_env) <- ~Long+Lat

####plot maps of the genetic summary stats on the landscape
#color palette for leaflet plot
pal_one <- colorNumeric(palette = "viridis", domain = NULL, na.color = "#00000000")


#plot genetic summary stats on a map and output the raster
for(stat in gen_sum_vec){
  #create raster of genetic summary stat per cell
  pi_raster_one <- rasterize(pi_env, sa_clim_1d, fun = "first", field = stat)
  writeRaster(pi_raster_one, filename = paste0("rasters_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".tif"), format = "GTiff")
  #leaflet plot. 
  l <- leaflet() %>% 
    addTiles() %>%
    addRasterImage(pi_raster_one, colors = pal_one,  opacity = 0.8) %>%
    addLegend(pal = pal_one, values = values(pi_raster_one))
  #write to plots folder
  mapview::mapshot(l, selfcontained = FALSE, url = paste0(getwd(), "/plots_", Sys.Date(), "/", str_replace_all(stat, "\\.", "_"), "_", Sys.Date(), ".html"))
}
 

#' 
#' 
