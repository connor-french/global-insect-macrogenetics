###Plotting bigass phylogeography data for RCN meeting
library(data.table)
library(raster)
library(tidyverse)
library(tidylog)
library(ggspatial)
library(rnaturalearth)

###########################################################

########## Species and individual counts ##################

###########################################################


#read in data
test_nuc <- fread("test_nuc_2019-05-21.csv") %>% 
  dplyr::select(-nucleotides, -markercode) %>% 
  rename(species = bin_uri, cell = cells) %>% 
  mutate(species = str_remove_all(species, "BOLD:"))


#make sure that my filtering worked and a maximum of 100 individuals per species per cell are included
test_nuc %>% 
  group_by(cell) %>% 
  count(species, sort = TRUE) %>% 
  top_n(20, n) %>% 
  View()

#also check for minimum of 3 individuals per species per cell
test_nuc %>% 
  group_by(cell) %>% 
  count(species, sort = TRUE) %>% 
  top_n(-1, n) %>% 
  View()

#check the number of alignments 
sp_counts <- test_nuc %>% 
  group_by(species) %>%
  count(cell, sort = TRUE) %>% 
  ungroup() %>% 
  count(cell, sort = TRUE) %>% 
  mutate(cell = as.factor(cell)) 


sp_counts %>%
  transmute(total = sum(n)) %>% 
  distinct(total)

#cell with the most species
sp_counts %>% 
  top_n(20, n) %>%
  mutate(cell = fct_reorder(cell, n)) %>% 
  ggplot(aes(x = cell, y = n)) +
  geom_col() +
  labs(x = "Cell number", y = "N") +
  coord_flip() +
  theme_minimal()

#cell with the least species. Cool, my filter to leave out any cell with less than 10 species worked
sp_counts %>% 
  top_n(-20, n) %>%
  mutate(cell = fct_reorder(cell, n)) %>% 
  ggplot(aes(x = cell, y = n)) +
  geom_col() +
  labs(x = "Cell number", y = "N") +
  coord_flip() +
  theme_minimal()

#ggplot of right-skewed distribution that looks like pi for illustrative purposes
pi_dist <- enframe(rbeta(100, 1, 100), value = "pi", name = NULL)
pi_plot <- pi_dist %>% 
  ggplot() +
  geom_histogram(aes(x = pi, y = ..density..), fill = "gold") +
  labs(y = "Frequency", x = "Avg. pairwise genetic distance per species") +
  theme_minimal() +
  theme(text = element_text(size=35, color = "white"),
        panel.background = element_rect(fill = "transparent", color = "white"), # bg of the panel
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

ggsave(pi_plot, 
       filename = "/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/Conferences/Evolution-2019/poster/pi_plot.png", 
       bg= "transparent",
       device = "png")



#plot rasters 
sp_rast <- raster("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/BigAss-bird-phylogeography/BigAss-phylogeography/rasters_2019-06-15/hill_one_avg_2019-06-15.tif")
sp_poly <- rasterToPolygons(sp_rast) %>% 
  sf::st_as_sf()

env_rast <- raster("/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/BigAss-bird-phylogeography/ndvi_dhi_combined-v6/ndvi_dhi_combo_downscaled_1-degree.tif")

world_map <- ne_coastline(returnclass = "sf")
mad_map <- ne_countries(country = "madagascar", returnclass = "sf")


mad_rast <- raster::crop(sp_rast, mad_map)
mad_poly <- rasterToPolygons(mad_rast) %>% 
  sf::st_as_sf()

mad_env <- crop(env_rast, mad_map)
mad_env_poly <- rasterToPolygons(mad_env) %>% 
  sf::st_as_sf()

sample_map <- ggplot() + 
  geom_sf(data = mad_env_poly, aes(fill = ndvi_dhi_combo_downscaled_1.degree)) +
  scale_fill_gradient(low = "white",high = "black") +
  labs(fill = "DHI", color = "white") +
  layer_spatial(mad_map, fill = "transparent", color = "white") +
  #ggnewscale::new_scale_fill() +
  #layer_spatial(sp_rast) +
  #scale_fill_viridis_c(na.value = "transparent", direction = -1) +
  #labs(fill = "Hill One", color = "white") +
  coord_sf(datum = NA) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        legend.text=element_text(color = "white"),
        legend.title = element_text(color = "white"),
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent", color = "white"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg)
        )
ggsave(sample_map, 
       filename = "/Users/connorfrench/Dropbox/Old_Mac/School_Stuff/CUNY/Conferences/Evolution-2019/poster/mad_env.png", 
       bg= "transparent",
       device = "png")

#################################
#########filter out extreme value with 12,000 species
################################
sp_rast_nomax <- sp_rast
max_index <- sp_rast_nomax %>% raster::which.max() 

#assign NA to huge value
values(sp_rast_nomax)[max_index] <- NA

#plot
ggplot() + 
  layer_spatial(world_map) +
  layer_spatial(sp_rast_nomax) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal()
  

##################################################
##########filter out cell with 1000 or more species
##################################################

#assign NA to values > 1000
sp_rast_mask_1000 <- sp_rast > 1000

sp_rast_1000 <- mask(sp_rast, sp_rast_mask_1000, maskvalue = TRUE)

#plot
ggplot() + 
  layer_spatial(world_map) +
  layer_spatial(sp_rast_1000) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal()

##################################################
##########filter out cell with 100 or more species
##################################################
#assign NA to values > 1000
sp_rast_mask_100 <- sp_rast > 100

sp_rast_100 <- mask(sp_rast, sp_rast_mask_100, maskvalue = TRUE)

#plot
ggplot() + 
  layer_spatial(world_map) +
  layer_spatial(sp_rast_100) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal()




###########################################################

################ Genetic diversity ########################

###########################################################

#read in sequence statistics dataframe. replace path if you're conducting this analysis on a different day than the first two steps
pi_df_file <- paste0("pi_df_one_2019-05-26.csv")
pi_df_one <- fread(pi_df_file) %>% 
  as_tibble() %>% 
  left_join(test_nuc, by = c("species", "cell")) %>% 
  distinct(species, cell, .keep_all = TRUE)

pi_df_one %>% 
  filter(avg_pi < 0.01) %>% 
  count(cell, sort = TRUE) %>% 
  filter(n > 9)

#filter for pi less than 0.01
pi_filter_hund <- pi_df_one %>% 
  filter(avg_pi < 0.01) %>% 
  group_by(cell) %>% 
  filter(n() > 9) %>% 
  ungroup()

#histogram of avg pi
pi_filter_hund %>% 
  ggplot(aes(x = avg_pi)) +
  geom_histogram(bins = 100)

#histogram of sd pi
pi_filter_hund %>% 
  ggplot(aes(x = sd_pi)) +
  geom_histogram(bins = 100)

#scatterplot of avg pi ~ sd pi
pi_filter_hund %>% 
  ggplot(aes(x = avg_pi, y = sd_pi)) +
  geom_point(alpha = 0.01) +
  geom_smooth() +
  geom_hline(yintercept = 0.01, col = "red") +
  theme_minimal()


pi_filter_hund %>% 
  filter(avg_pi > 0) %>% 
  count()








