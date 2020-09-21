###############################################################################
## Project:     VAST covariates across exptraolation grid
## Author:      Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description: Assign bathymetry value from the EFH data lyaer to each grid in
##              the Gulf of Alaska extrapolation grid.
###############################################################################
rm(list = ls())

##################################################
####   Import packages
##################################################
library(rgdal)
library(raster)

##################################################
####   Set up directories
##################################################
github_dir <- "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/data/"

##################################################
####   Import EFH bathymety raster and extrapolation grid
####   Import CPUE data
##################################################
bathy <- raster::raster(
  paste0(github_dir, "EFH_bathymetry/aigoa_bathp1c/dblbnd.adf"))
grid = read.csv("extrapolation_grid/GOAThorsonGrid.csv")
data = read.csv(paste0(github_dir, "GOA_multspp.csv"))



##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = sp::SpatialPointsDataFrame(
  coords = grid[, c("Longitude", "Latitude")],
  data = data.frame(area = grid$Shape_Area),
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
grid_shape_aea = sp::spTransform(x = grid_shape,
                                 CRSobj = crs(bathy))
grid_shape_aea@data$depth =  raster::extract(x = bathy,
                                             y = grid_shape_aea,
                                             method = "simple")

##################################################
####   There are 3048 gridd cells without depths, For these stations, we 
####   extract raster values averaged within a 15 km buffer. This seemed to 
####   be the lowest radius found to provide a value for each unknown station
##################################################
distance = 200
how_many_NAs = sum(is.na(grid_shape_aea@data$depth))
while (how_many_NAs) {
  grid_shape_aea@data$depth[is.na(grid_shape_aea@data$depth)] <- 
    raster::extract(x = bathy,
                    y = grid_shape_aea[is.na(grid_shape_aea@data$depth),],
                    buffer = distance,
                    na.rm = T,
                    fun = mean)
  
  how_many_NAs = sum(is.na(grid_shape_aea@data$depth))
  print(distance)
  print(how_many_NAs)
  distance = distance + 200
}

##################################################
####   Plot depth covariate of the extrapolation grid
##################################################
spplot(grid_shape_aea[, "depth"], 
       col.regions = rev(terrain.colors(1000)),
       pch = 16, 
       cex = 0.1,
       cuts = 100,
       colorkey = T)

Extrapolation_depths = grid
Extrapolation_depths$Depth = grid_shape_aea@data$depth

Extrapolation_depths[, c("E_km", "N_km")] <- project(
  xy = coordinates(Extrapolation_depths[, c("Longitude", "Latitude")]), 
  proj = "+proj=utm +zone=5N" )

##################################################
####   Save
##################################################
save(list = "Extrapolation_depths", 
     file = paste0(github_dir, 'Extrapolation_depths.RData'))
