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
library(furrr)
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

print("Read in data")

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

print("Read in climate data")

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
#' 
#' Write the sequence data in the BOLD csv to nexus files. Writing nexus files for species where there are at least three observations and occupy cells where there are more than nine species per cell.  Each nexus file is labeled as "speciesname.cell.nex".
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

#density plot of number of species per cell
count_density_sp <- test_nuc %>% 
  group_by(cells) %>% 
  count() %>% 
  ggplot(aes(x = n, y = ..density..)) +
  geom_density(fill = "grey", alpha = 0.4) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7) +
  labs(x = "Number of individuals", y = "Density") +
  theme_minimal() 

ggsave(filename = paste0("plots_", Sys.Date(), "/filtered_species_counts_per_cell_", Sys.Date(), ".pdf"), count_density_sp)

#convert to a spatial data frame
coordinates(test_nuc) <- ~Long+Lat 

#make a raster out of the counts of species per cell
count_pts_1d_sp <- rasterize(test_nuc, sa_clim_1d, fun = function(x, ...) {length(unique(x))}, field = "bin_uri") 

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

print("Filtered sequence data")

#write to output folder
fwrite(test_nuc, file = paste0("output_", Sys.Date(), "/test_nuc_", Sys.Date(), ".csv"))
#test_nuc <- fread("bold-seqs-insects-10.txt")


#split data frame by cell number and species
test_nuc <- as.data.frame(test_nuc)
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

#split into sublists because the computer gets gassed
files <- split(files, cut(seq_along(files), 5, labels = FALSE))

print("Listed nexus files")

#only run if I need to re-calculate pi for everything. Takes forever
plan(multicore, workers = 8)

for (file in files){
  #calculate pi statistics for each species within each cell
  pi_df_one <- file %>% furrr::future_map_dfr(gen_dist_calc) 
  ###Write this to a csv. Can read in later if you don't want to run all of the stats again.
  fwrite(pi_df_one, file = paste0("output_", Sys.Date(), "/pi_df_one_", stringr::str_replace(Sys.time(), " ", "_"), ".csv"))
}

print("Calculated genetic summary statistics")

