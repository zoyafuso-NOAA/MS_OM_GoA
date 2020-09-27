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
bathy = bathy + abs(min(values(bathy), na.rm = T))

goa_grid = read.csv(paste0(github_dir, "extrapolation_grid/GOAThorsonGrid.csv"))
goa_grid <- goa_grid[, c("Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2 
names(goa_grid) <- c("Area_km2", "Lon", "Lat")

data = read.csv(paste0(github_dir, "GOA_multspp.csv"))

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = sp::SpatialPointsDataFrame(
  coords = goa_grid[, c("Lon", "Lat")],
  data = goa_grid,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
grid_shape_aea = sp::spTransform(x = grid_shape,
                                 CRSobj = crs(bathy))
grid_shape_aea@data$depth =  raster::extract(x = bathy,
                                             y = grid_shape_aea,
                                             method = "simple")

##################################################
####   Plot depth covariate of the extrapolation grid
##################################################
spplot(grid_shape_aea[, "depth"], 
       col.regions = rev(terrain.colors(1000)),
       pch = 16, 
       cex = 0.1,
       cuts = 100,
       colorkey = T)

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
Extrapolation_depths <- goa_grid
Extrapolation_depths[, c("E_km", "N_km")] <- project(
  xy = coordinates(Extrapolation_depths[, c("Lon", "Lat")]), 
  proj = "+proj=utm +zone=5N +units=km" )
Extrapolation_depths$DEPTH_EFH = grid_shape_aea@data$depth

Extrapolation_depths_subset = subset(x = Extrapolation_depths,
                                     subset = DEPTH_EFH < max(data$DEPTH_EFH) &
                                       DEPTH_EFH > min(data$DEPTH_EFH) )

vast_grid = FishStatsUtils::gulf_of_alaska_grid
vast_grid[, c("E_km", "N_km")] <- project(
  xy = coordinates(vast_grid[, c("Lon", "Lat")]), 
  proj = "+proj=utm +zone=5N +units=km" )

{png(paste0(github_dir, "Subsetted_Extrapolation_Grid.png"),
    width = 6,
    height = 5,
    units = "in",
    res = 500)

par(mfrow = c(2, 1), mar = c(0,0,0,0))
plot(N_km ~ E_km,
     data = Extrapolation_depths,
     pch = 15,
     cex = 0.25,
     col = 'red', 
     axes = F)
points(N_km ~ E_km,
       data = Extrapolation_depths_subset,
       col = 'black',
       cex = 0.25,
       pch = 15)
mtext(side = 1, line = -2, text = "Proposed New Grid")

plot(N_km ~ E_km,
     data = Extrapolation_depths,
     pch = 15,
     cex = 0.25,
     col = 'red',
     axes = F)
points(N_km ~ E_km,
       data = vast_grid,
       col = 'black',
       cex = 0.25,
       pch = 15)
mtext(side = 1, line = -2, text = "FishStatsUtils Grid")

dev.off()}

##################################################
####   scale grid bathymetry values to standard normal, using the mean and sd
####   of the BTS data
##################################################
BTS_mean <- mean(log(data$DEPTH_EFH))
BTS_sd   <- sd(log(data$DEPTH_EFH))

Extrapolation_depths <- Extrapolation_depths_subset
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
