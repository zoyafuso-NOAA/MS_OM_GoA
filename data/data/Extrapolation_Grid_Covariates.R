##################################
## Assign bathymetry value to each extrapolation grid
## Gulf of Alaska
#################################
setwd("/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/data/")

##################################
## Import Libraries
#################################
library(marmap); library(sp); library(RANN); library(raster);

##################################
## Load Extrapolation Grid used in VAST
#################################
load('Spatial_Settings.RData')

##################################
## Extract fine-scale bathymetry map 
## Convert latlon to UTM zone 5
#################################
xmin <- -170
xmax <- -132
ymin <- 52
ymax <- 60.5
bathymap <- getNOAA.bathy(lon1 = xmin, lon2 = xmax,
                          lat1 = ymin, lat2 = ymax,
                          resolution = 1)
bathymap <- fortify.bathy(bathymap)

bathymap_coord = sp::SpatialPoints(coords = bathymap[,c('x', 'y')],
                                   proj4string = CRS('+proj=longlat') )
cord.UTM <- sp::spTransform(bathymap_coord, CRS("+proj=utm +zone=5N"))
bathymap[,c('E_km', 'N_km')] = cord.UTM@coords / 1000

#####################
## Assign bathymetry values for each extrapolation grid cell to
## the nearest point in the bathymetry map
#####################

bathy_idx = RANN::nn2(query=Extrapolation_List$Data_Extrap[,c('E_km','N_km')],
                      data = bathymap[,c('E_km', 'N_km')],
                      k = 300)$nn.idx

Extrapolation_List$Data_Extrap$depth = -bathymap$z[bathy_idx[,1]]

#####################
## Plot locations where depths are negative (land?) 
#####################
plot(Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], pch = '.')
points(Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$depth <=0,c('E_km', 'N_km')], pch = 16, col = 'red')

#####################
## Assign nearest bathymetry values to those grid cells where depths are 
## negative. This involes an iterative while loop
#####################

neg_depths = sum(Extrapolation_List$Data_Extrap$depth <=0)
k = 2
while(neg_depths != 0){
  idxs = which(Extrapolation_List$Data_Extrap$depth <= 0)
  Extrapolation_List$Data_Extrap$depth[idxs] = -bathymap$z[bathy_idx[idxs,k]]
  neg_depths = sum(Extrapolation_List$Data_Extrap$depth <=0)
  k = k + 1
}

#############################
## Center depth and calculate depth^2
#############################
Extrapolation_List$Data_Extrap$DEPTH = scale(Extrapolation_List$Data_Extrap$depth)

#############################
## Plot Bathyetry Field
#############################
test = raster::rasterize(
  x = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')],
  y = raster(nrows=100, ncols=320,
             xmn=min(Extrapolation_List$Data_Extrap$E_km),
             xmx=max(Extrapolation_List$Data_Extrap$E_km),
             ymn=min(Extrapolation_List$Data_Extrap$N_km),
             ymx=max(Extrapolation_List$Data_Extrap$N_km),
             crs = CRS("+proj=utm +zone=5N")),
  field = Extrapolation_List$Data_Extrap$depth)
plot(test, axes = F)

##########################
## Save
##########################
