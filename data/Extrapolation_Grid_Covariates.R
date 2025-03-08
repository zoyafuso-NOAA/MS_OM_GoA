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
library(sp)
library(rgeos)

##################################################
####   Set up directories
##################################################
EFH_dir <- "G:/Oyafuso/data/aigoa_bathp1c/"
github_dir <- "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/"
github_dir2 <- "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/data/" 

##################################################
####   Import EFH bathymety raster
####   Import Current Strata (need to match the current spatial domain)
####   Import Extrapolation grid
####   Import CPUE data
####   Untrawlable areas
##################################################
bathy <- raster::raster( paste0(EFH_dir, "dblbnd.adf"))
bathy <- bathy + abs(min(values(bathy), na.rm = T))

current_survey_strata <- rgdal::readOGR(
  paste0(github_dir, "shapefiles/goa_strata.shp"))

current_survey_mask <- rgdal::readOGR(
  paste0(github_dir, "shapefiles/goagrid_polygon.shp"))
current_survey_mask <- sp::spTransform(x = current_survey_mask,
                                       CRSobj = crs(bathy))
current_survey_mask <- rgeos::gUnaryUnion(spgeom = current_survey_mask)

goa_grid <- read.csv(paste0(github_dir, 
                            "extrapolation_grid/GOAThorsonGrid.csv"))

goa_grid <- goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2 
names(goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

data = read.csv(paste0(github_dir, "GOA_multspp.csv"))

goa_grid_nountrawl <- read.csv(
  paste0(github_dir, 
         "extrapolation_grid/GOA_ALL_nountrawl.csv"))

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = sp::SpatialPointsDataFrame(
  coords = goa_grid[, c("Lon", "Lat")],
  data = goa_grid,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

grid_shape_aea = sp::spTransform(x = grid_shape,
                                 CRSobj = crs(bathy))
grid_shape_aea@data$DEPTH_EFH =  raster::extract(x = bathy,
                                                 y = grid_shape_aea,
                                                 method = "simple")

##################################################
####   Remove cells not already in the goa stratification
##################################################
grid_shape_aea <- raster::intersect(x = grid_shape_aea,
                                    y = current_survey_mask)

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
grid_shape_aea <- subset(grid_shape_aea,
                         DEPTH_EFH >= min(data$DEPTH_EFH) & 
                           DEPTH_EFH <= max(data$DEPTH_EFH))

##################################################
####   Plot depth covariate of the extrapolation grid
##################################################
spplot(grid_shape_aea[, "DEPTH_EFH"], 
       col.regions = rev(terrain.colors(1000)),
       pch = 16, 
       cex = 0.1,
       cuts = 100,
       colorkey = T)

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
Extrapolation_depths <- grid_shape_aea@data
Extrapolation_depths[, c("E_km", "N_km")] <- project(
  xy = coordinates(Extrapolation_depths[, c("Lon", "Lat")]), 
  proj = "+proj=utm +zone=5N +units=km" )

##################################################
####   Create indices to easily subset <700 m cells and untrawlable cells
##################################################
Extrapolation_depths$shallower_than_700m <- cells_shallower_than_700m <-
  Extrapolation_depths$DEPTH_EFH <= 700
Extrapolation_depths$trawlable <- cells_trawlable <-
  Extrapolation_depths$Id %in% goa_grid_nountrawl$Id
Extrapolation_depths$shallow_trawlable <- 
  rowSums(Extrapolation_depths[, c("shallower_than_700m", 
                                   "trawlable")]) == 2

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
####   Add current strata labels to each grid cell
##################################################
Extrapolation_depths$stratum <- raster::extract( x = current_survey_strata, 
                                                 y = grid_shape_aea)$STRATUM
Extrapolation_depths$stratum[is.na(Extrapolation_depths$stratum)] <- 0

##################################################
####   Save
##################################################
save(list = c("Extrapolation_depths"),
     file = paste0(github_dir, "Extrapolation_depths.RData"))

save(list = c("Extrapolation_depths"),
     file = paste0(github_dir2, "Extrapolation_depths.RData"))