library(raster)
library(here)
library(tictoc)
library(tidyverse)

dir.create("~/data")
dir.create("~/data/climate_agg")


out_path <- here("data", "climate_agg")
res_list <- list(high = 96.5, medium = 193, low = 385.9)

crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

# transform world raster to high resolution
template_high <- raster(here("data", "template_high.tif"))


# transform world raster to medium resolution
template_medium <- raster(here("data", "template_medium.tif"))

# function to resample rasters to any equal area grid of the three resolutions we're considering.


# crs: define the crs you want to project to. It must be an equal area projection. The default crs is behrmann equal area. 
# x: a raster or raster stack
# km: km resolution you want (area will be km x km)
# is_categorical: if is_categorical is true, using "ngb" interpolation
resample_equal_area <- function(x, 
                                crs = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs", 
                                km,
                                is_categorical = FALSE) {
  
  # project raster to new coordinate system (this has no values- we put values back in later)
  x_ext <- raster::projectExtent(x, crs)
  
  # make the resolution square since we're interested in an equal area projection
  raster::res(x_ext) <- raster::xres(x_ext)
  
  # add original values to our projection. We can't take shortcuts like doing this instead of doing projectExtent first because it removes cells at the extreme latitudes.
  x_repro <- raster::projectRaster(x, x_ext)
  
  if (km == 96.5) {
    end_raster = template_high
  } else if (km == 193) {
    end_raster = template_medium
  } else if (km == 385.9) {
    end_raster = template_low
  } else stop("Incorrect resolution specification. Must be either 96.5 or 193")
  
  # resample our raster of interest to the beginning raster  
  if (is_categorical) {
    m <- "ngb"
  } else m <- "bilinear"
  
  
  x_resampled <- raster::resample(x_repro, end_raster, method = m) 
  
  
  return(x_resampled)
}


# function to write aggregated rasters to file, using a specified prefix for the file and layer names
write_raster <- function (x, prefix) {
  names(x) <- paste0(prefix, "_", names(x))
  out_file <- file.path(out_path, paste0(names(x), ".tif"))
  writeRaster(x, filename = out_file, overwrite = TRUE)
}




soil_path <- dir(here("data"), pattern = "*.tif$", full.names = TRUE)
soil_rast_list <- purrr::map(soil_path, raster)

tic()
soil_agg_high <- purrr::map(soil_rast_list, ~resample_equal_area(.x, km = res_list$high, is_categorical = TRUE))
toc()
# write rasters to file
invisible(purrr::map(soil_agg_high, ~write_raster(.x, prefix = "soil_high")))

tic()
soil_agg_medium <- purrr::map(soil_rast_list, ~resample_equal_area(.x, km = res_list$medium, is_categorical = TRUE))
toc()
# write rasters to file
invisible(purrr::map(soil_agg_medium, ~write_raster(.x, prefix = "soil_medium")))





