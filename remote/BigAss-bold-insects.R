#' ---
#' title: "BigAss BOLD"
#' always_allow_html: yes
#' output: # use rmarkdown::render(<your-rmd-file.rmd>, output_format ="all" in the console to render both outputs
#'   html_notebook:
#'     theme: flatly
#'     highlight: tango
#'   github_document:
#'     toc: true
#' ---
#' Load packages, read in data, and filter for NAs
## ----include = FALSE-----------------------------------------------------
library(data.table)
library(tidyverse)
library(tidylog)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)
#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
source("R/pad_short_seqs_function.R") #pad short sequences with Ns
source("R/nexus_write_function.R") #write nexus files 
source("R/align_calc_gen_dist.R") #align seqs and calculate mean pi per species per cell
source("R/hill_num_calc.R") #calculate hill number of pi. Note* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.
source("R/pi_summary_function.R")
source("R/sumstat_plot_function.R")

#make output folder
dir.create(paste0("output_", Sys.Date()))
#make plots folder
dir.create(paste0("plots_", Sys.Date()))
#make raster folder
dir.create(paste0("rasters_", Sys.Date()))



#read in data and remove any row that contains NA values for latitude, longitude, and bin. Only reading in the first 20000 rows for debugging purposes
bold.data <- fread("../bold_data_insects.txt", na.strings = c("","NA"),  quote = "", verbose = TRUE) %>% 
  filter_at(vars(lat, lon, bin_uri), all_vars(!is.na(.))) 
  


#convert to sf object for plotting
coord_points <- st_as_sf(bold.data, coords = c("lon", "lat"), 
                         crs = 4326, agr = "constant")

#' 
#' 
#' Make a raster of the number of individuals per cell we've obtained after filtering for at least three individuals per species per cell and plot. The outlier individuals get filtered out at the next filtering step.
## ------------------------------------------------------------------------

#Reading in an environmental raster (world clim 2.0 bio1 raster) at 10 arc-minute resolution to use as a raster that utilizes lat-long coordinates for its cells. The raster isn't what's important, what is important is the raster resolution. 
f <- list.files("../wc2.0_10m_bio/", full.names = TRUE) #raster location


sa_clim_1d <- stack(f) %>% #read in bioclim rasters, downscale the resolution to 1 degree, and crop them
#  crop(bounds) %>% 
  aggregate(fact = 6) 
  

#get the cell number of each coordinate pair for filtering.
pts_ext_1d <- raster::extract(sa_clim_1d, coord_points, fun = "count", sp = TRUE, cellnumbers = TRUE) %>% 
  as.data.frame() %>%
  setnames(old = c("coords.x1", "coords.x2"), new = c("Long", "Lat")) %>% #rename coordinates 
  mutate(cells = as.factor(cells)) %>% #need cells numbers as factors so I can count their frequency
  group_by(bin_uri, cells) %>% #group the data set by species, then by cell number
  filter(n() > 2) %>% #retain only observations where there are more than two species observations per cell
  ungroup() %>%
  drop.levels()

#density plot of number of species per cell
count_density <- pts_ext_1d %>% 
  group_by(cells) %>% 
  count() %>% 
  ggplot(aes(x = n, y = ..density..)) +
  geom_density(fill = "grey", alpha = 0.4) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7) +
  labs(x = "Number of individuals", y = "Density") +
  theme_minimal() 

ggsave(filename = paste0("plots_", Sys.Date(), "/unfiltered_species_counts_per_cell_", Sys.Date(), ".pdf"), count_density)

coordinates(pts_ext_1d) <- ~Long+Lat #convert to a spatial data frame


count_pts_1d <- rasterize(pts_ext_1d, sa_clim_1d, fun = "count", field = "bin_uri") #make a raster out of the counts of observations per cell

pal.vert <- colorBin(palette = "inferno", bins = 10, domain = NULL, pretty = TRUE, na.color = "#00000000")


#plot count map
all_plot <- leaflet(data = sa_clim_1d) %>% 
  addTiles() %>%
  addRasterImage(count_pts_1d, colors = pal.vert,  opacity = 0.8) %>%
  addLegend(pal = pal.vert, values = values(count_pts_1d))

#write raster to file
writeRaster(count_pts_1d, filename = paste0("rasters_", Sys.Date(), "/all_species_count_", Sys.Date(),".tif"), format = "GTiff")

mapview::mapshot(all_plot, selfcontained = FALSE, url = paste0(getwd(), "/plots_", Sys.Date(), "/all_species_count_", Sys.Date(),".html"))

#' 
#' Plot of the number of species per cell after filtering for cells that contain ten or more species
## ------------------------------------------------------------------------
#convert to data frame for filtering
pts_ext_1d_sp <- as.data.frame(pts_ext_1d) %>% #convert to data frame for filtering
  group_by(cells) %>%
  filter(n_distinct(bin_uri) > 9) %>%
  ungroup()

#write to output file
fwrite(pts_ext_1d_sp, file = paste0("output_", Sys.Date(), "/pts_ext_1d_sp_", Sys.Date(), ".csv"))

#density plot of number of species per cell
count_density_sp <- pts_ext_1d_sp %>% 
  group_by(cells) %>% 
  count() %>% 
  ggplot(aes(x = n, y = ..density..)) +
  geom_density(fill = "grey", alpha = 0.4) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7) +
  labs(x = "Number of individuals", y = "Density") +
  theme_minimal() 

ggsave(filename = paste0("plots_", Sys.Date(), "/filtered_species_counts_per_cell_", Sys.Date(), ".pdf"), count_density_sp)

#convert to a spatial data frame
coordinates(pts_ext_1d_sp) <- ~Long+Lat 

#make a raster out of the counts of species per cell
count_pts_1d_sp <- rasterize(pts_ext_1d_sp, sa_clim_1d, fun = function(x, ...) {length(unique(x))}, field = "bin_uri") 

pal.species <- colorBin(palette = "inferno", bins = 10, domain = NULL, pretty = TRUE, na.color = "#00000000")

#plot count map
filter_count_map <- leaflet(data = sa_clim_1d) %>% 
  addTiles() %>%
  addRasterImage(count_pts_1d_sp, colors = pal.species,  opacity = 0.8) %>%
  addLegend(pal = pal.species, values = values(count_pts_1d_sp))

#write raster to file
writeRaster(count_pts_1d_sp, filename = paste0("rasters_", Sys.Date(), "/filter_species_count_", Sys.Date(),".tif"), format = "GTiff")

#write to plots folder
mapview::mapshot(filter_count_map, selfcontained = FALSE, url = paste0(getwd(), "/plots_", Sys.Date(), "/filter_species_count_", Sys.Date(), ".html"))

#' 
#' 
#' Write the sequence data in the BOLD csv to nexus files. Writing nexus files for species where there are at least three observations and occupy cells where there are at least nine species per cell.  Each nexus file is labeled as "speciesname.cell.nex".
## ------------------------------------------------------------------------
#####filter data for the ideal number of individuals and species per cell
test_nuc <- pts_ext_1d %>% 
  raster::as.data.frame(xy = TRUE) %>%
  dplyr::select(recordID, bin_uri, cells, markercode, nucleotides, Lat, Long, contains("bio_10m")) %>% 
  distinct(recordID, .keep_all = TRUE) %>% #only retain unique individuals
  na.omit() %>%
  group_by(bin_uri, cells) %>% #group the data set by species, then by cell number
  filter(str_detect(markercode, "COI"), !str_detect(markercode, "COII")) %>% #Filter for only COI sequences
  filter(n() > 2) %>% #retain only cells where there are more than two species observations per cell
  ungroup() %>%
  group_by(cells) %>%
  filter(n_distinct(bin_uri) > 9) %>% #retain cells with 10 or more species
  ungroup() %>%
  drop.levels()

#write to output folder
fwrite(test_nuc, file = paste0("output_", Sys.Date(), "/test_nuc_", Sys.Date(), ".csv"))
#test_nuc <- fread("bold-seqs-insects-10.txt")


#split data frame by cell number and species
species_seq_split_one <- test_nuc %>%
  as.data.frame() %>% #some functions don't like spatial data frames
  drop.levels() %>%
  split(.$cells) %>% #split into a list of data frames, grouped by cell. Splitting by both cells and species at once doesn't work.
  lapply(drop.levels) %>% #drop any levels in the data frame. Have to perform first because extra factor levels can mess up the split function
  lapply(function(x){
    split(x, x$bin_uri)}) %>% #split each cell by species
  drop.levels()

#create directory to put nexus files
nexus_folder <- paste0("output_", Sys.Date(), "/bold_seqs_", Sys.Date())
dir.create(nexus_folder)

#write nexus files to output folder
nexus_write_fun(species_seq_split_one, out_folder = nexus_folder)

print("Wrote nexus files")

#' 
#' 
#' Calculate genetic diversity statistics for each cell. I had to manually edit several nexus files due to some sequences denoting gaps with spaces and some using dashes.
## ------------------------------------------------------------------------
#get a list of the nexus files
files <- list.files(nexus_folder, full.names = TRUE)

print("Listed nexus files")

#only run if I need to re-calculate pi for everything. Takes forever
pi_df_one <- files %>% purrr::map_dfr(gen_dist_calc) #calculate avg pi and sd pi for each species within each cell

print("Calculated genetic summary statistics")

###Write this to a csv. Can read in later if you don't want to run all of the stats again.
fwrite(pi_df_one, file = paste0("output_", Sys.Date(), "/pi_df_one_", Sys.Date(), ".csv"))
#pi_df_one <- fread("raw-pi-insects-10.csv")

print("Wrote summary stats to file")

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
