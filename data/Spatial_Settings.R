###############################
## Spatial Settings for VAST
###############################
rm(list = ls())

library(VAST)

setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/')

## Import Data
data = read.csv(file = 'data/GOA_multspp.csv')
data = subset(x = data, YEAR > 1990)

# Prepare the Data-frame for catch-rate data
Data_Geostat = data.frame( "spp"=data$SCI, 
                           "Year"=data$YEAR, 
                           "Catch_KG"=data$BIOMASS, 
                           "AreaSwept_km2"=data$EFFORT, 
                           "Vessel"=0, 
                           "Lat"=data$LATITUDE, 
                           "Lon"=data$LONGITUDE,
                           "DEPTH" = data$DEPTH,
                           "DEPTH2" = data$DEPTH2)

#Filter Species
Data_Geostat = subset(x = Data_Geostat, 
                      spp %in% c('Gadus macrocephalus', 'Atherestes stomias',
                                 'Sebastes alutus', 'Gadus chalcogrammus',
                                 'Hippoglossus stenolepis', 'Solea solea',
                                 'Glyptocephalus zachirus'))

## Spatial settings: The following settings define the spatial resolution 
## for the model, and whether to use a grid or mesh approximation
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
#grid_size_km = 50
n_x = 250   # Specify number of stations (a.k.a. "knots")

## Stratification for results
strata.limits <- data.frame(
  'STRATA' = c("All_areas", "west_of_140W"),
  'west_border' = c(-Inf, -Inf),
  'east_border' = c(Inf, -140)
)

## Extrapolation grid: We also generate the extrapolation grid appropriate 
## for a given region.

Extrapolation_List = make_extrapolation_info( Region= "Gulf_of_Alaska", 
                                              strata.limits = strata.limits )

## Derived objects for spatio-temporal estimation: And we finally generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`

Spatial_List = make_spatial_info( grid_size_km=sqrt(13.72), 
                                  n_x=n_x, 
                                  Method=Method, 
                                  Lon=Data_Geostat[,'Lon'], 
                                  Lat=Data_Geostat[,'Lat'], 
                                  Extrapolation_List=Extrapolation_List,
                                  DirPath=getwd(), 
                                  Save_Results=TRUE )

# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, 
                      "knot_i"=Spatial_List$knot_i )

#plot data points and knots
plot(Spatial_List$latlon_i[,2:1])
points(Spatial_List$latlon_x[,2:1], col = 'red', pch = 16, cex = 1)

save.image('Spatial_Settings.RData')
