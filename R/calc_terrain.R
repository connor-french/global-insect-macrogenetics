library(raster)

gdem <- raster(here("data", "climate_raw", "terrain", "elevation.tif"))

# calculate raster values
terrain <- terrain(gdem, opt = c("slope", 
                                 "aspect",
                                 "TPI",
                                 "TRI",
                                 "roughness"), 
                   unit = "degrees")

filenames <- here::here("data", "climate_raw", "terrain", names(terrain))

writeRaster(terrain, filenames, bylayer = TRUE, format = "GTiff")
