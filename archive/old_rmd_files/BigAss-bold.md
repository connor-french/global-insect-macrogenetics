BigAss BOLD
================

Load packages, read in data, and filter for NAs

Make a raster of the number of individuals per cell we've obtained after filtering for at least three individuals per species per cell and plot. The outlier individuals get filtered out at the next filtering step.

``` r
#Reading in an environmental raster (world clim 2.0 bio1 raster) at 10 arc-minute resolution to use as a raster that utilizes lat-long coordinates for its cells. The raster isn't what's important, what is important is the raster resolution. 
f <- list.files("/Users/connorfrench/Dropbox/Old_Mac/climate-data/wc2.0_10m_bio/", full.names = TRUE) #raster location

bounds <- extent(-85, -30, -45, 15) #Reduce extent of raster to our area of interest- South America

sa_clim_1d <- stack(f) %>% #read in bioclim rasters, downscale the resolution to 1 degree, and crop them
  aggregate(fact = 6) %>%
  crop(bounds)



#get the cell number of each coordinate pair for filtering.
pts_ext_1d <- raster::extract(sa_clim_1d, coord_points, fun = "count", sp = TRUE, cellnumbers = TRUE) %>% 
  as.data.frame() %>%
  setnames(old = c("coords.x1", "coords.x2"), new = c("Long", "Lat")) %>% #rename coordinates 
  mutate(cells = as.factor(cells)) %>% #need cells numbers as factors so I can count their frequency
  group_by(species_name, cells) %>% #group the data set by species, then by cell number
  filter(n() > 2) %>% #retain only observations where there are more than two species observations per cell
  ungroup() %>%
  drop.levels()

coordinates(pts_ext_1d) <- ~Long+Lat #convert to a spatial data frame

count_pts_1d <- rasterize(pts_ext_1d, sa_clim_1d, fun = "count", field = "species_name") #make a raster out of the counts of observations per cell

pal.vert <- colorBin(palette = "inferno", bins = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,300,500,1000,1500,2000), domain = NULL, pretty = TRUE, na.color = "#00000000")


#plot count map
leaflet(data = sa_clim_1d) %>% 
  addTiles() %>%
  addRasterImage(count_pts_1d, colors = pal.vert,  opacity = 0.8) %>%
  addLegend(pal = pal.vert, values = values(count_pts_1d))
```

    ## Warning in colors(.): Some values were outside the color scale and will be
    ## treated as NA

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-4ddf5de2ce295fe5b274">{"x":{"options":{"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}}},"calls":[{"method":"addTiles","args":["//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",null,null,{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false,"attribution":"&copy; <a href=\"http://openstreetmap.org\">OpenStreetMap<\/a> contributors, <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\">CC-BY-SA<\/a>"}]},{"method":"addRasterImage","args":["data:image/png;base64,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",[[15,-85],[-45,-30]],0.8,null,null,null]},{"method":"addLegend","args":[{"colors":["#000004","#0D0829","#280B54","#480B6A","#65156E","#82206C","#9F2A63","#BB3655","#D44842","#E8602D","#F57D15","#FB9E07","#FAC127","#F3E55C","#FCFFA4"],"labels":["0 &ndash; 10","10 &ndash; 20","20 &ndash; 30","30 &ndash; 40","40 &ndash; 50","50 &ndash; 60","60 &ndash; 70","70 &ndash; 80","80 &ndash; 90","90 &ndash; 100","100 &ndash; 300","300 &ndash; 500","500 &ndash; 1,000","1,000 &ndash; 1,500","1,500 &ndash; 2,000"],"na_color":null,"na_label":"NA","opacity":0.5,"position":"topright","type":"bin","title":null,"extra":null,"layerId":null,"className":"info legend","group":null}]}],"limits":{"lat":[-45,15],"lng":[-85,-30]}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
Plot of the number of species per cell after filtering for cells that contain ten or more species

``` r
pts_ext_1d_sp <- as.data.frame(pts_ext_1d) %>% #convert to data frame for filtering
  group_by(cells) %>%
  filter(n_distinct(species_name) > 9) %>%
  ungroup()


coordinates(pts_ext_1d_sp) <- ~Long+Lat #convert to a spatial data frame

count_pts_1d_sp <- rasterize(pts_ext_1d_sp, sa_clim_1d, fun = function(x, ...) {length(unique(x))}, field = "species_name") #make a raster out of the counts of species per cell

pal.species <- colorBin(palette = "inferno", bins = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,200,300,500, 1000, 2000), domain = NULL, pretty = TRUE, na.color = "#00000000")

#plot count map
leaflet(data = sa_clim_1d) %>% 
  addTiles() %>%
  addRasterImage(count_pts_1d_sp, colors = pal.species,  opacity = 0.8) %>%
  addLegend(pal = pal.species, values = values(count_pts_1d_sp))
```

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-e1fb4bf25072e1fd93ca">{"x":{"options":{"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}}},"calls":[{"method":"addTiles","args":["//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",null,null,{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false,"attribution":"&copy; <a href=\"http://openstreetmap.org\">OpenStreetMap<\/a> contributors, <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\">CC-BY-SA<\/a>"}]},{"method":"addRasterImage","args":["data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADcAAAA8CAYAAADCHCKFAAACgUlEQVRoge2av07bQBzHfyGpGuWU0QygSOmEU8QTIAY63QNYnVErLxbyG7DAG1CpS6V269SNhQ0hRF+grWBEQtAIj5Ej+occ00XmsH127s7nS+6zRLLj831//+7uJwNYLJZa4yKPsNe67cGza8Zg9ORVkubpqlji/aHbHpD9/i7J8x7PsxgFxHfCykVyxQEAfI8msNpcz/TC6P6ikXbdRR7pwzL0Oi24Hv9/ZgTfCQlGgTLRqZNKQie02lwHAIDL+Bv3mSQYBaTXaUEPTeBrNIQ+LE/vXcEd/PhwBC/e/Sk1ZlFaeTep5c8fToRecja+BRgD0LFG9xcNjAKy1VlRJowLRgERLQgu8qah120PhMcTJi03RMekAmshDqOAxHuO0etbZrwnhaGDSE9ezBPJ/FTGv88vcxdtlbA5KTVH4z1Haq4VMZJMAdz3fXn9vnKv0d2Kjm2aNuyJxJKD9p2FKuZSWDKZFyaxVR8ytbPf39UqbqHWQyW4yJN+AK1VrleVf9qFinqRFcB6sYzAQq29sohYmG0TZrUNiyBdXNnWXxGM2ETwul/U41llvhZhmUfSs74TTsXSX4wCcja+nd73nZBQUZvN7VJhKlUcr2qyE/sUHTYAnna1r+AObh5+wU/8hrx99fvJ88fxR32NqlmbuGy4ddsDs07k7GT/nq6l5o/M5lQlOVfUCzSsRcp/5Qx3NshwZyNzcaZgFGhrKQrBtgvTRJQt9UbgIo+oEFb5OpeHUbm2MOg+uReibG4YISrJXFY1EXQYo9JqWet2ILW+KSFZ2HMyhNXeKCYVk1rtUGRjxZnKTOI2m9v1Lw7A+WqPxXdC0kMTuI6XYBTZHbzFMgOPyQ8ZmDG0eYMAAAAASUVORK5CYII=",[[15,-85],[-45,-30]],0.8,null,null,null]},{"method":"addLegend","args":[{"colors":["#000004","#0D0829","#280B54","#480B6A","#65156E","#82206C","#9F2A63","#BB3655","#D44842","#E8602D","#F57D15","#FB9E07","#FAC127","#F3E55C","#FCFFA4"],"labels":["0 &ndash; 10","10 &ndash; 20","20 &ndash; 30","30 &ndash; 40","40 &ndash; 50","50 &ndash; 60","60 &ndash; 70","70 &ndash; 80","80 &ndash; 90","90 &ndash; 100","100 &ndash; 200","200 &ndash; 300","300 &ndash; 500","500 &ndash; 1,000","1,000 &ndash; 2,000"],"na_color":null,"na_label":"NA","opacity":0.5,"position":"topright","type":"bin","title":null,"extra":null,"layerId":null,"className":"info legend","group":null}]}],"limits":{"lat":[-45,15],"lng":[-85,-30]}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
Count the number of phyla per cell to characterize the higher-order taxonomic diversity per cell. This could influence calculations of genetic diversity.

``` r
count_pts_1_ph <- rasterize(pts_ext_1d_sp, sa_clim_1d, fun = function(x, ...) {length(unique(x))}, field = "phylum_name") #make a raster out of the counts of species per cell
pal.phylum <- colorBin(palette = "inferno", bins = c(1,2,3,4,5,6,7,8,9,10), domain = NULL, pretty = TRUE, na.color = "#00000000")


#plot phylum map
leaflet(data = sa_clim_1d) %>% 
  addTiles() %>%
  addRasterImage(count_pts_1_ph, colors = pal.phylum,  opacity = 0.8) %>%
  addLegend(pal = pal.phylum, values = values(count_pts_1_ph))
```

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-a694d0665c8d184f45df">{"x":{"options":{"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}}},"calls":[{"method":"addTiles","args":["//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",null,null,{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false,"attribution":"&copy; <a href=\"http://openstreetmap.org\">OpenStreetMap<\/a> contributors, <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\">CC-BY-SA<\/a>"}]},{"method":"addRasterImage","args":["data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADcAAAA8CAYAAADCHCKFAAACMElEQVRoge2aMU7DMBSG/yKPGSIkmHGpeoAOuQEdsqOIvQPq0JUz9AJcgr1DuUEHZgYgA0MXxoyVzBCZJmnSOOmL7QZ/U5Wm1vv9/J6dPwUcDofVcC8UKtfOhrMOvluYsYm5qLuBe6GI/IU4nr06Aazm991QKw4ANrsPAEDkLyoC3A3Kr6eiuTdFnKxRnAQ5carBNkVJHAAEbNRi+FR0KiwlzSATUuhsvG0xrhpHxWWX0uRy95fBU0iF7jN9t3mpyHrHcC8kqJV9lvKfDVEURNMMsgL1ULoksmLiZGVm2RBQGXhfBFqEhpqkaSRtOdwLyYamFqa2SdNlqzZ2M1ljwoqtQi//SqyjGT1+YO3lus+K6qXAMnrfss2K63EP0EQX59DseMYzpC8A43V+WgCHE1Xs2EYFdhlAs7GVrT11qjzM9iyHcwEAkT9H5M+ph6fiuPMsv6u6p975zqPZG2Eim1kZaJysBsWg5bWAjTAbb/H2cwUAePp6Vo6ZeFnW1UN+yUrjKfIXQrrScbJGwEZYDufi8fom9+smwjqgXTNJrQiW29u6fIdATpnRW1Y/3w8Tsm7bQbc8RLUJcC8U75+3oOq4WsRVUWb2ysZBgRZxphxrS2xyJrg3BUA7EUaXZRH3TsJGjD+XqdF0ozb+DNYM7oVCntQdAF6Dezv/h0KF1XW3D+486kg5czTCrJ8U4yaNMladUKhx4rJwbwp5yLWdEzJnf901Epffp+j9SYfjH/ALoL3cAiZH+LkAAAAASUVORK5CYII=",[[15,-85],[-45,-30]],0.8,null,null,null]},{"method":"addLegend","args":[{"colors":["#000004","#210C4A","#57106E","#89226A","#BB3655","#E35932","#F98D0A","#F9C932","#FCFFA4"],"labels":["1 &ndash; 2","2 &ndash; 3","3 &ndash; 4","4 &ndash; 5","5 &ndash; 6","6 &ndash; 7","7 &ndash; 8","8 &ndash; 9","9 &ndash; 10"],"na_color":null,"na_label":"NA","opacity":0.5,"position":"topright","type":"bin","title":null,"extra":null,"layerId":null,"className":"info legend","group":null}]}],"limits":{"lat":[-45,15],"lng":[-85,-30]}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
Write the sequence data in the BOLD csv to nexus files. Writing nexus files for species where there are at least three observations and occupy cells where there are at least nine species per cell. Each nexus file is labeled as "speciesname.cell.nex".

``` r
#####filter data for the ideal number of individuals and species per cell
test.nuc <- as.data.frame(pts_ext_1d) %>%
  dplyr::select(recordID, species_name, cells, marker_codes, nucleotides, contains("bio_10m")) %>% 
  distinct(recordID, .keep_all = TRUE) %>% #only retain unique individuals
  na.omit() %>%
  group_by(species_name, cells) %>% #group the data set by species, then by cell number
  filter(str_detect(marker_codes, "COI"), !str_detect(marker_codes, "COII")) %>% #Filter for only COI sequences
  filter(n() > 2) %>% #retain only cells where there are more than two species observations per cell
  ungroup() %>%
  group_by(cells) %>%
  filter(n_distinct(species_name) > 9) %>% #retain cells with 10 or more species
  ungroup() %>%
  drop.levels()



#####Code for padding short sequences with gaps
#function to pad the ends of shorter sequences with "N"s so all sequences are the same length for alignments
gap.fix.fun <- function(nuc) {
  nuc <- nuc[order(sapply(nuc, length), decreasing = TRUE)] #order sequences in descending order by length. This removes the necessity to do a bunch of unnecessary pairwise comparisons among sequences. Hashing this out because I'm ordering the data frame by sequence length before I apply this function. 
  nuc.out <- vector("list", length(nuc)) #establish an empty vector to fill
  nuc.out[[1]] <- nuc[[1]] #replace the first sequence in the out vector
  names(nuc.out) <- names(nuc) #keep the names of the input vector
  for (i in seq(2, length(nuc)))
    if(length(nuc[[i]]) <= length(nuc[[1]])) #if the sequence is shorter than the longest sequence, append the missing data designator "N" to the end of the sequence
      nuc.out[[i]] <- append(nuc[[i]], rep(c("N"), length(nuc[[1]]) - length(nuc[[i]])))
  return(nuc.out)
  }




#split data frame by cell number and species
species.seq.split.one <- test.nuc %>%
  as.data.frame() %>% #some functions don't like spatial data frames
  drop.levels() %>%
  split(.$cells) %>% #split into a list of data frames, grouped by cell. Splitting by both cells and species at once doesn't work.
  lapply(drop.levels) %>% #drop any levels in the data frame. Have to perform first because extra factor levels can mess up the split function
  lapply(function(x){
    split(x, x$species_name)}) %>% #split each cell by species
  drop.levels()

#create directory to put nexus files
dir.create("../bold-seqs-10")
```

    ## Warning in dir.create("../bold-seqs-10"): '../bold-seqs-10' already exists

``` r
#function to write nexus files
nexus.write.fun <- function(x) {
  for (i in seq_along(x)) {
    df <- x[[i]]
    for (j in seq_along(df)) {
      nuc.df <- as.data.frame(df[[j]]) #loop through each species per cell
      #nuc.df <- nuc.df[order(sapply(nuc.df[,5], length), decreasing = TRUE),] #order df by nucleotide length (need to do this for padding short sequences with "N"s)
      nuc.vec <- strsplit(nuc.df[,5], "") #nucleotide column. Need to split nucleotides into individual characters.
      names(nuc.vec) <- paste(nuc.df[,1]) #names are the recordIDs
      nuc.pad <- gap.fix.fun(nuc.vec) #pad nucleotides with "N"s
      write.nexus.data(nuc.pad, file = paste0("../bold-seqs-10/", str_replace(unique(nuc.df[,2]), " ", "-"), ".", unique(nuc.df[,3]), ".nex"), missing = "N", interleaved = FALSE) #write to a nexus file, which is named by the species and cell number
    }
  }
}

nexus.write.fun(species.seq.split.one)
```

A few functions to calculate mean pi per species per cell, the number of sequences per species per cell, and the second Hill number per cell. I'm aligning sequences with clustal omega after reading them in. Note\* using clustal omega in the ape R package requires having a copy of clustal omega downloaded and accessible either through their PATH or have it indicated in the clustalomega() argument.

``` r
#create function to read, align, and calculate mean raw genetic distance (pi) from the sequences.
gen.calc.fun <- function(x) {
  #print(x) #print if debugging
  rawseq <- read.nexus.data(x)
  binseq <- as.DNAbin(rawseq)
  alignseq <- clustalomega(binseq) 
  dist <- dist.dna(alignseq, model = "raw")
  return(mean(dist))
}


## Get one hill number from a list of genetic distances. Original python code written by Isaac Overcast
hill.calc <- function(dists, order) { 
  if (order == 0) {
    return(length(dists[dists > 0]))
  }
  if (order == 1) {
    h1 = exp(entropy::entropy(dists))
    return(h1)
  }
  else {
    tot = sum(dists)
    proportions = dists[dists > 0]/tot
    prop_order = proportions**order
    h2 = sum(prop_order)**(1/(1-order))
    return(h2)
  }
}
```

Calculate genetic diversity statistics for each cell. I had to manually edit several nexus files due to some sequences denoting gaps with spaces and some using dashes.

``` r
#get a list of the nexus files
files <- paste0("../bold-seqs-10/", list.files("../bold-seqs-10"))

#only run if I need to re-calculate pi for everything. Takes forever
pi.total.one <- files %>% sapply(gen.calc.fun) #calculate pi for all each species within each cell

#function to only extract cell number from the name of the sequence. Some species names include a number, so I need to extract only the last number grouping
sec.num.fun <- function(n) {
  for (i in seq_along(n)) {
    if (length(n[[i]]) > 1){
      n[[i]] <- n[[i]][length(n[[i]])]
    }
  }
  return(n)
}

pi.names.one <- str_match_all(names(pi.total.one), "[0-9]+") %>% sec.num.fun() %>% unlist() #make a vector of the cell numbers for each pi calculation

pi.df.one <- bind_cols(cells = pi.names.one, pi = unname(pi.total.one)) #create dataframe of pi calculations
###I wrote this to a csv so I don't have to re-run the pi calculation script
#write.csv(pi.df.one, file = "raw-pi-10.csv")
#pi.df.one <- read.csv("raw-pi-10.csv")


#function to calculate pi for each cell number. I'm using dplyr::summarise because plyr also has a summarise function that can create a lot of confusion.
pi.summary.fun <- function(df) {
  gdf <- group_by(df, cells)
  sdf <- sample_n(gdf, 10)
  sum.df <- dplyr::summarise(sdf, median.pi = median(pi), mean.pi = mean(pi), sd.pi = sd(pi), hill.zero = hill.calc(pi, 0), hill.one = hill.calc(pi, 1), hill.two = hill.calc(pi, 2), shannon = entropy::entropy(pi))
  return(sum.df)
}


#replicate the sampling 1000 times
sum.list<- replicate(1000, pi.summary.fun(pi.df.one), simplify = FALSE) 

#summarise these samples into a final df
sum.df <- bind_rows(sum.list) %>%
  group_by(cells) %>% 
  dplyr::summarise(median.pi.avg = mean(median.pi), median.pi.sd = sd(median.pi), mean.pi.avg = mean(mean.pi),  mean.pi.sd = sd(mean.pi), hill.zero.avg = mean(hill.zero), hill.zero.sd = sd(hill.zero), hill.one.avg = mean(hill.one), hill.one.sd = sd(hill.one), hill.two.avg = mean(hill.two), hill.two.sd = sd(hill.two), shannon.avg = mean(shannon), shannon.sd = sd(shannon))
write.csv(sum.df, file = "pi-summary-10.csv")

#make new data frame including the pi values.
pi.plot <- merge(pts_ext_1d_sp, sum.df, by = "cells") %>% 
  raster::as.data.frame(xy = TRUE) %>% #have to convert to data frame to omit NAs
  subset(select = c(recordID, species_name, cells, marker_codes, median.pi.avg, median.pi.sd, mean.pi.avg, mean.pi.sd, hill.zero.avg, hill.zero.sd, hill.one.avg, hill.one.sd, hill.two.avg, hill.two.sd, shannon.avg, shannon.sd, Long, Lat)) %>%
  na.omit() %>%
  drop.levels()

#convert back to a spatial data frame
coordinates(pi.plot) <- ~Long+Lat


#create new raster of mean pi values per cell
pi.raster.one <- rasterize(pi.plot, sa_clim_1d, fun = "first", field = "mean.pi.avg") 

#color palette for leaflet plot
pal.one <- colorNumeric(palette = "viridis", domain = NULL, na.color = "#00000000")

#leaflet plot
leaflet() %>% 
  addTiles() %>%
  addRasterImage(pi.raster.one, colors = pal.one,  opacity = 0.8) %>%
  addLegend(pal = pal.one, values = values(pi.raster.one))
```

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-0125ffea8cbb58eba604">{"x":{"options":{"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}}},"calls":[{"method":"addTiles","args":["//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",null,null,{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false,"attribution":"&copy; <a href=\"http://openstreetmap.org\">OpenStreetMap<\/a> contributors, <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\">CC-BY-SA<\/a>"}]},{"method":"addRasterImage","args":["data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADcAAAA8CAYAAADCHCKFAAADUElEQVRogWNgGAWjYBQMauAhlfN/oN0wCugBHJ3aByymmYhRpJ/dR5EDzWJ6/1sHdtPdk0R5joGBgcEksfd/ztkokhxo69f1n+n3X4Y/HIwMv7mZMGLRQ6mYph5mJFahSWLv/zPzi4lWDwPOti3/XxtxMTD8Y2DgffqX4RcPE8MfLkYG0RMfGP6zMzO81eNlODe7iGRziQEEDY05kfx/icVcsi2XX9jxX0nmFcPnxdIMLD//MzD+YWA4vrqEJp5BB3iTpUp73/9HLeoUWfD/DyPD660yDGcWFDP+iX3LcGTSTIrMoxh4SOei5g2xDIrzhk4RpFCy2lVGt4IFa/LwkM79/8lSnuEPJyPDqcWk57PBArAmyx1PJzOyf/jNwMDAwKDWSFk1MAqgQLG/97/c7C7aBqZ2Sd9/d/6kAYkx3QLUlPLlmRz13KE4ofe/Yn8vXT3mLpRKNfvcuGLhZmEtLFTXNP2/HVJH14LEJgDRPDuyoXTIFmKjYBRQC1hE9AzP+o6apdmgAW7siP6bO2/C8PMgNmAZ2vPfJIG+dSFdgb17x4B6zo07bvgGLl2AbmHff7Um6vYc1OsR5tG8IU0IuFi30MUBRimkBSLRo1/4wJ6jNfC2oJ1nJ0Ue9RDPRNHvbIsacPpbagYuJu08O/+jd1uGFXC2a6Wq5xz3FQ6NwPLQqCDo0ILz4VjVaJX3/dcqGyKpwtWiCcWhenl9/5HrT7OY3v8TrzvD+aSmCKoUKDDgIZJGkuW7T+DuELsZ1f1n+f6PIV9zL1zN3kPVA9eJ9RBJ++/CGEJWspGfgTpR4mZUNzSSHwMDA4OHbD6KY239uv6vuG2C4QFn25b/HsolVPEYVZMlLuAunPr//9dvBNXp5UMKix13e6iS/FioYQhRAM25z6xYGCJUz6CIXppI3dkeusTczrezGRm5uehh1SigK7DzGeCW/yjAA9y1q4Zf7BgnDeOxllEwHAB6L3tQAfmp9F8FRAkguoVimNH3X+oApHXkZE9eTxt9lcSgA2axvf9N44dGCUhyw/mHECPDf7q0SCkHJHuO6TcDA8OQiDcyegWcr/8xMA5XzzEwMDAI3vzBoJc3+EehSE6WXC9/MjD8/89waRJtlhGOglEwzAEAC2wQCsHqdSYAAAAASUVORK5CYII=",[[15,-85],[-45,-30]],0.8,null,null,null]},{"method":"addLegend","args":[{"colors":["#440154 , #443A83 16.8750674823219%, #2C728E 37.6895496039318%, #20A486 58.5040317255418%, #76D054 79.3185138471517%, #FDE725 "],"labels":["0.005","0.010","0.015","0.020"],"na_color":null,"na_label":"NA","opacity":0.5,"position":"topright","type":"numeric","title":null,"extra":{"p_1":0.168750674823219,"p_n":0.793185138471517},"layerId":null,"className":"info legend","group":null}]}],"limits":{"lat":[-45,15],"lng":[-85,-30]}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
NEED TO FIGURE OUT CELL NUMBERING PROBLEM
