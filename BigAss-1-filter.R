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
#library(furrr)
library(leaflet)
library(sf)
library(raster)
library(ape)
library(gdata)
library(entropy)
library(ggspatial)
library(rnaturalearth)

#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
source("R/pad_short_seqs_function.R") #pad short sequences with Ns
source("R/fasta_write_function.R") #write fasta files 
source("R/calc_gen_dist.R") #align seqs and calculate mean pi per species per cell
source("R/hill_num_calc.R") #calculate hill number of pi. Note* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.
source("R/pi_summary_function.R")
source("R/sumstat_plot_function.R")
source("R/equal_area_resampling.R")

#### Set theme for the plots in this project ----
theme_set(theme_minimal())


#### Make output folders ----

#if you want to work out of a folder from an earlier date, replace this string with the date
todays_date <- Sys.Date()

#folder for the entire project's output to go into
todays_results <- paste0("results_", todays_date)
dir.create(todays_results)

#output folder for spreadsheets, tables, etc.
output <- paste0(todays_results, "/output")
dir.create(output)

#folder for plots
plots <- paste0(todays_results, "/plots")
#make plots folder
dir.create(plots)

#folder for rasters
rasters <- paste0(todays_results, "/rasters")
#make raster folder
dir.create(rasters)


###Start filtering
#read in data and remove any row that contains NA values for latitude, longitude, and bin. 
bold_data <- fread("../bold_data_insects.txt", na.strings = c("","NA"),  quote = "") %>% 
  filter_at(vars(lat, lon, bin_uri), all_vars(!is.na(.))) 

print("Read in data")

#convert to sf object for plotting and project to behrmann equal area crs
behrmann_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
coord_points <- st_as_sf(bold_data, coords = c("lon", "lat"), 
                         crs = 4326, agr = "constant") %>% 
  st_transform(crs = behrmann_crs)

#' 
#' 
#' Make a raster of the number of individuals per cell we've obtained after filtering for at least three individuals per species per cell and plot. The outlier individuals get filtered out at the next filtering step.
## ------------------------------------------------------------------------

###read in environmental raster. I have to hard code this because I can't host them on the github and they need to be downloaded
f <- raster("/Users/connorfrench/Dropbox/Old_Mac/climate-data/chelsa_10min/10min/bio_1.tif")

#resample the env't to 100 km cells. 
sa_clim_1d <- f %>% 
  resample_equal_area(km = 100)

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

### Num of individual plots, unfiltered dataset ----
#histogram of number of individuals per cell
count_histogram <- pts_ext_1d %>% 
  as.data.frame(xy = TRUE) %>% 
  group_by(cells) %>% 
  count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7, bins = 50) +
  labs(x = "Number of individuals", y = "Frequency") 

ggsave(filename = paste0(todays_results,"/plots/unfiltered_species_counts_per_cell.pdf"), 
       count_histogram)

#convert to a spatial data frame
coordinates(pts_ext_1d) <- ~Long+Lat 


#make a raster out of the counts of observations per cell
count_pts_1d <- rasterize(pts_ext_1d,
                          sa_clim_1d,
                          fun = "count",
                          field = "bin_uri")


##plot individual count map
#first grab world base map
world_map <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
  st_transform(crs = behrmann_crs)

#make a ggplot of the map
all_map <- ggplot() +
  geom_sf(data = world_map) +
  ggspatial::layer_spatial(data = count_pts_1d) +
  scale_fill_viridis_c(na.value = NA) +
  labs(title = "Unfiltered individuals counts", 
       subtitle = paste0("# cells = ", 
                         ncell(na.omit(values(count_pts_1d))), 
                         "; # individuals = ", 
                         sum(na.omit(values(count_pts_1d))),
                         "; # species = ",
                         pts_ext_1d %>% as_tibble() %>% filter(!duplicated(bin_uri)) %>% count()), 
       fill = "Individuals count") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


#write raster to file
writeRaster(
  count_pts_1d,
  filename = paste0(todays_results, "/rasters/all_species_count.tif"),
  format = "GTiff"
)

#save plot to file
ggsave(
  plot = all_map,
  filename = paste0(todays_results, "/plots/all_species_count.pdf")
)
  


#' 
#' 
#' Write the sequence data in the BOLD csv to fasta files. Writing fasta files for species where there are at least three observations and occupy cells where there are more than nine species per cell.  Each fasta file is labeled as "speciesname.cell.nex".
## test_nuc ----
#####filter data for the ideal number of individuals and species per cell.
test_nuc <- pts_ext_1d %>%
  raster::as.data.frame(xy = TRUE) %>%
  distinct(recordID, .keep_all = TRUE) %>% #only retain unique individuals
  filter_at(vars(recordID, bin_uri, cells, markercode, nucleotides, Lat, Long),
            all_vars(!is.na(.))) %>%
  group_by(bin_uri, cells) %>% #group the data set by species, then by cell number
  filter(str_detect(markercode, "COI"),
         !str_detect(markercode, "COII")) %>% #Filter for only COI sequences
  filter(n() > 2) %>% #retain only cells where there are more than two species observations per cell
  sample_n(if (n() < 100)
    n()
    else
      100) %>% #subsample individuals if there are 100 or more individuals per cell. This keeps sample sizes within an order of magnitude
  ungroup() %>%
  group_by(cells) %>%
  filter(n_distinct(bin_uri) > 9) %>% #retain cells with 10 or more species
  ungroup() %>%
  drop.levels()

###Num of individuals per cell, filtered data set ----
###histogram of number of individuals per cell
count_histogram_filtered <- test_nuc %>% 
  group_by(cells) %>% 
  count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7, bins = 50) +
  labs(x = "Number of individuals", y = "Frequency") 

ggsave(filename = paste0(todays_results, "/plots/filtered_individual_counts_per_cell.pdf"), count_histogram_filtered)

###barplot of the number of individuals from each order in the dataset
order_barplot <- test_nuc %>% 
  count(order_name) %>% 
  mutate(order_name = fct_reorder(order_name, n)) %>% 
  ggplot(aes(x = order_name, y = n)) +
  geom_col() +
  labs(x = "Order", y = "# of species") +
  coord_flip()

ggsave(
  plot = order_barplot,
  filename = paste0(todays_results, "/plots/order_individuals_barplot.pdf")
)


#Num of species per cell, filtered data set ----
###histogram of number of individuals per cell
count_histogram_filtered_sp <- test_nuc %>% 
  group_by(cells) %>% 
  filter(!duplicated(bin_uri)) %>%
  count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram(fill = "darkgreen", color = "black", alpha = 0.7, bins = 50) +
  labs(x = "Number of species", y = "Frequency") 

ggsave(filename = paste0(todays_results, "/plots/filtered_species_counts_per_cell.pdf"), count_histogram_filtered_sp)




###barplot of the number of species from each order in the dataset
order_barplot_sp <- test_nuc %>% 
  group_by(bin_uri) %>% 
  filter(!duplicated(order_name)) %>% 
  count(order_name, sort = TRUE) %>% 
  ungroup() %>% 
  count(order_name) %>% 
  mutate(order_name = fct_reorder(order_name, n)) %>% 
  ggplot(aes(x = order_name, y = n)) +
  geom_col() +
  labs(x = "Order", y = "# of species") +
  coord_flip()

ggsave(
  plot = order_barplot_sp,
  filename = paste0(todays_results, "/plots/order_species_barplot.pdf")
)

#convert to a spatial data frame
coordinates(test_nuc) <- ~Long+Lat 

#make a raster out of the counts of species per cell
count_pts_1d_sp <- rasterize(test_nuc, sa_clim_1d, fun = function(x, ...) {length(unique(x))}, field = "bin_uri") 


###make a ggplot of the map
filter_count_map_sp <- ggplot() +
  geom_sf(data = world_map) +
  ggspatial::layer_spatial(data = count_pts_1d_sp) +
  scale_fill_viridis_c(na.value = NA) +
  labs(title = "Filtered species counts", 
       subtitle = paste0("# cells = ", ncell(na.omit(values(count_pts_1d_sp))), "; # species = ", sum(na.omit(unique(values(count_pts_1d_sp))))), 
       fill = "Species count") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#write raster to file
writeRaster(count_pts_1d_sp, filename = paste0(todays_results, "/rasters/filter_species_count.tif"), format = "GTiff")

ggsave(
  plot = filter_count_map_sp,
  filename = paste0(todays_results, "/plots/filtered_species_map.pdf")
)

###barplot of the top 20 counts of number of individuals per species per cell. I'm trying to see if there are any species that have an insane number of individuals per cell for a particular cell.
count_ind <- test_nuc %>% 
  as_tibble() %>% 
  group_by(bin_uri, cells) %>% 
  count(sort = TRUE) %>% 
  mutate(bin_cells = str_c(bin_uri, cells, sep = "_")) %>%
  mutate(bin_cells = str_remove_all(bin_cells, "BOLD:")) %>% 
  ungroup() %>% 
  top_n(20, n) %>% 
  mutate(bin_cells = fct_reorder(bin_cells, n)) %>% 
  ggplot(aes(x = bin_cells, y = n)) +
  geom_col() +
  labs(x = "Species bin_Cell number", y = "N") +
  coord_flip() +
  theme_minimal()

ggsave(filename = paste0(todays_results, "/plots/top_20_sp_counts.pdf"), count_ind)


print("Filtered sequence data")

#write to output folder
test_nuc <- as.data.frame(test_nuc)
fwrite(test_nuc, file = paste0(todays_results, "/output/test_nuc.csv"))


