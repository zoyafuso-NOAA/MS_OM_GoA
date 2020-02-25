###############################
## Spatial Settings for VAST
###############################
rm(list = ls())

library(VAST)
setwd( 'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

modelno = '4b'
if(!dir.exists(paste0(getwd(), '/VAST_output', modelno, '/'))) {
  dir.create(paste0(getwd(), '/VAST_output', modelno, '/'))
}

## Import Data
data = read.csv(file = 'data/data/GOA_multspp.csv')

# Prepare the Data-frame for catch-rate data
Data_Geostat = data.frame( "spp"=data$SPECIES_NAME,
                           "Year"=data$YEAR,
                           "Catch_KG"=data$WEIGHT,
                           "AreaSwept_km2"=data$EFFORT,
                           "Vessel"=0,
                           "Lat"=data$LATITUDE,
                           "Lon"=data$LONGITUDE,
                           "DEPTH" = data$DEPTH,
                           "DEPTH2" = data$DEPTH2)
rm(data)

#Drop factor levels of unused Species
Data_Geostat$spp = droplevels(Data_Geostat$spp)

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

fine_scale = T
Spatial_List = make_spatial_info( n_x=n_x, 
                                  Method=Method, 
                                  Lon_i=Data_Geostat[,'Lon'],
                                  Lat_i=Data_Geostat[,'Lat'],
                                  Extrapolation_List=Extrapolation_List,
                                  DirPath=getwd(), 
                                  fine_scale = fine_scale,
                                  Save_Results=T )

# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, 
                      "knot_i"=Spatial_List$knot_i )

#plot data points and knots
plot(Spatial_List$latlon_i[,2:1])
points(Spatial_List$latlon_x[,2:1], col = 'red', pch = 16, cex = 1)

save.image(paste0(getwd(), '/VAST_output', modelno, '/Spatial_Settings.RData') )
