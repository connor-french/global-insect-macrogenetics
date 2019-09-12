### I typically run this on a cluster since it takes forever otherwise. On 11 cores, takes ~20 minutes

# this is the script I'm using to resample the values to an equal area projection
source("R/equal_area_resampling.R")

### read in environmental data for plotting
# We have three different data sets: 19 CHELSA bioclims (chelsa), 3 DHI variables (dhi), and global habitat heterogeneity (ghh) from Tuanmu and Jetz 2015
chelsa_dir <- "~/Desktop/connorfrench/bigass-phylogeography/10min"
dhi_file <- "~/Desktop/connorfrench/bigass-phylogeography/ndvi_dhi_combined-v6/ndvi_dhi_combo_small.tif" #only one file, so don't need to list contents
#ghh_dir <- "~/Desktop/connorfrench/bigass-phylogeography/global_habitat_heterogeneity"

chelsa_files <- list.files(chelsa_dir, pattern = ".tif$", full.names = TRUE)
#ghh_files <- list.files(ghh_dir, pattern = ".tif$", full.names = TRUE)



chelsa <- stack(chelsa_files)
dhi <- stack(dhi_file)
#ghh <- stack(ghh_files)

#parallelizing the process and recording the time because I'm curious
start_time <- Sys.time()
raster::beginCluster()

### aggregate to equal area projection at chosen resolution in km
res <- 200

chelsa_agg <- resample_equal_area(chelsa, km = res)
writeRaster(chelsa_agg, filename = "~/Desktop/connorfrench/bigass-phylogeography/chelsa_300km.tif", format = "GTiff")
dhi_agg <- resample_equal_area(dhi, km = res)
writeRaster(dhi_agg, filename = "~/Desktop/connorfrench/bigass-phylogeography/dhi_300km.tif", format = "GTiff")
#ghh_agg <- resample(ghh, km = res)
#writeRaster(ghh_agg, filename = "~/Desktop/connorfrench/bigass-phylogeography/ghh_100km/tif", format = "GTiff")

raster::endCluster() 
end_time <- Sys.time()
end_time - start_time

### combine into a single raster stack
#clim <- stack(chelsa_agg, dhi_agg)

#writeRaster(clim, 
#            filename = "~/Desktop/connorfrench/bigass-phylogeography/all_clim_200km.tif",
#            format = "GTiff")
