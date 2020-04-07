source("~/Desktop/connorfrench/bigass-phylogeography/equal_area_resampling.R")

# list all raster directories in master raster directory
raster_dirs <- list.dirs("~/Desktop/connorfrench/bigass-phylogeography/rasters/", full.names = TRUE)

# vector of raster resolutions
res_vec <- c(100, 200, 300)

# loop through directories and resolutions, writing out each resolution
for (dir in raster_dirs) {
  for (res in res_vec) {
    tifs <- list.files(dir, pattern = ".tif$", full.names = TRUE)
    r_stack <- raster::stack(tifs)
    
    # use the cluster to speed things up
    raster::beginCluster()
    
    # resample the rasters
    r_agg <- resample_equal_area(r_stack, km = res)
    raster::writeRaster(r_agg, file = paste0(dir, "/agg_", res, ".tif"))
    raster::endCluster() 
  }
  
}
