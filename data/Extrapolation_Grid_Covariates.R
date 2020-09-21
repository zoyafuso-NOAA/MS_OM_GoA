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
library(FishStatsUtils)

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
# grid = read.csv("extrapolation_grid/GOAThorsonGrid.csv")
grid = FishStatsUtils::gulf_of_alaska_grid
data = read.csv(paste0(github_dir, "GOA_multspp.csv"))

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = sp::SpatialPointsDataFrame(
  coords = grid[, c("Lon", "Lat")],
  data = grid,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
grid_shape_aea = sp::spTransform(x = grid_shape,
                                 CRSobj = crs(bathy))
grid_shape_aea@data$depth =  raster::extract(x = bathy,
                                             y = grid_shape_aea,
                                             method = "simple")

##################################################
####   There are 551 grid cells without depths, For these stations, we extract 
####   raster values averaged within an iteratively increasing buffer radius.  
####   The farthest radius was 13.4 km. 
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

Extrapolation_depths <- grid
Extrapolation_depths[, c("E_km", "N_km")] <- project(
  xy = coordinates(Extrapolation_depths[, c("Lon", "Lat")]), 
  proj = "+proj=utm +zone=5N" )
Extrapolation_depths$DEPTH_EFH = grid_shape_aea@data$depth

##################################################
####   Remove cells that are shallower and deeper than observed
####   Show in red which cells were removed
##################################################
Extrapolation_depths = subset(x = Extrapolation_depths,
                              subset = DEPTH_EFH < max(data$DEPTH_EFH) &
                                DEPTH_EFH > min(data$DEPTH_EFH) )
# Extrapolation_depths_subset = subset(x = Extrapolation_depths,
#                                      subset = DEPTH_EFH < max(data$DEPTH_EFH) &
#                                        DEPTH_EFH > min(data$DEPTH_EFH) )

# png(paste0(github_dir, "Subsetted_Extrapolation_Grid.png"),
#     width = 6,
#     height = 5,
#     units = "in",
#     res = 500)
# plot(N_km ~ E_km,
#      data = Extrapolation_depths,
#      pch = 15,
#      cex = 0.25,
#      col = 'red')
# points(N_km ~ E_km,
#        data = Extrapolation_depths_subset,
#        col = 'black',
#        cex = 0.25,
#        pch = 15)
# dev.off()

##################################################
####   scale grid bathymetry values to standard normal, using the mean and sd
####   of the BTS data
##################################################
BTS_mean <- mean(log(data$DEPTH_EFH))
BTS_sd   <- sd(log(data$DEPTH_EFH))

Extrapolation_depths$LOG_DEPTH_EFH <- log(Extrapolation_depths$DEPTH_EFH)
Extrapolation_depths$LOG_DEPTH_EFH_CEN <- 
  (Extrapolation_depths$LOG_DEPTH_EFH - BTS_mean) / BTS_sd
Extrapolation_depths$LOG_DEPTH_EFH_CEN_SQ <- 
  Extrapolation_depths$LOG_DEPTH_EFH_CEN ^ 2

##################################################
####   Save
##################################################
save(list = "Extrapolation_depths", 
     file = paste0(github_dir, 'Extrapolation_depths.RData'))
