#function to resample a raster to an equal area grid of a specified km resolution. currently limited to global-scale rasters.
#only tested from high to low resolutions, not the other way around.

#This function requires the raster package for raster manipulation and rnatural earth package for creating a blank world map
library(raster)
library(rnaturalearth)

#crs: define the crs you want to project to. It must be an equal area projection. The default crs is behrmann equal area. 
#x: a raster or raster stack
#km: km resolution you want (area will be km x km)
resample_equal_area <- function(x, crs = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs", km) {
  
  #project raster to new coordinate system (this has no values- we put values back in later)
  x_ext <- raster::projectExtent(x, crs)
  
  #make the resolution square since we're interested in an equal area projection
  raster::res(x_ext) <- raster::xres(x_ext)
  
  #add original values to our projection. We can't take shortcuts like doing this instead of doing projectExtent first because it removes cells at the extreme latitudes.
  x_repro <- raster::projectRaster(x, x_ext)
  
  #read in a world shape file. Might make this able to use countries as well
  world_shp <- rnaturalearth::ne_coastline(scale = "small")
  
  #convert the coordinate system to behrmann
  world_shp <- sp::spTransform(world_shp, crs)
  
  #turn it into a raster
  world_raster <- raster::raster(world_shp)
  
  #transform the raster to a certain resolution (in meters)
  raster::res(world_raster) <- km * 1000
  
  #resample our raster of interest to the world raster
  x_resampled <- raster::resample(x_repro, world_raster)
  
  return(x_resampled)
}


