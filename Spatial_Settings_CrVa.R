###############################
## Spatial Settings for VAST
###############################
rm(list = ls())

###################################
## Set up directories
###################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2, 'VM' = 3)[2]

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output6j')

## Import Data
data = read.csv(file = paste0(github_dir, 'data/data/GOA_multspp.csv') )

# Prepare the Data-frame for catch-rate data
Data_Geostat = data.frame( "spp"=data$SPECIES_NAME,
                           "Year"=data$YEAR,
                           "Catch_KG"=data$WEIGHT,
                           "AreaSwept_km2"=data$EFFORT,
                           "Vessel"=0,
                           "Lat"=data$LATITUDE,
                           "Lon"=data$LONGITUDE,
                           "LOG_DEPTH" = data$DEPTH, #centered log_depth
                           "LOG_DEPTH2" = data$DEPTH^2 )
rm(data)

#Drop factor levels of unused Species
spp_df = read.csv(paste0(github_dir, "spp_df.csv"), check.names=F, header = T, 
                  row.names = 'modelno')

which_spp = unlist(spp_df['6j',])

Data_Geostat = subset(Data_Geostat, spp %in% names(which_spp)[which_spp])
Data_Geostat$spp = droplevels(Data_Geostat$spp)

################################
## Assign 10 fold partitions of the data
################################
# Generate partitions in data
n_fold = 10

#Sort Data_Geostat
Data_Geostat = Data_Geostat[order(Data_Geostat$Year, Data_Geostat$spp),]
Data_Geostat$id = 1:nrow(Data_Geostat)

set.seed(2342)
foldno = lapply(X = split.data.frame(Data_Geostat, f = Data_Geostat$Year),
                FUN = function(test) {
                  row.idx = matrix(data = test$id, ncol = 15)
                  fold_no = sample(x = 1:n_fold, 
                                   size = nrow(row.idx), 
                                   replace = T)
                  return(split(row.idx, fold_no))
                })

test = lapply(X = foldno,
              FUN = function(test){
                lapply(test, FUN = function(x) Data_Geostat$id[x])
              })

for(iyear in names(test)){
  for(ifold in paste(1:10)){
    Data_Geostat[test[[iyear]][[ifold]],'fold'] = as.integer(ifold)
  }
}

## Spatial settings: The following settings define the spatial resolution 
## for the model, and whether to use a grid or mesh approximation
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
#grid_size_km = 50
n_x = 350   # Specify number of stations (a.k.a. "knots")

## Stratification for results
strata.limits <- data.frame(
  'STRATA' = c("All_areas"),#, "west_of_140W"),
  'west_border' = c(-Inf),#, -Inf),
  'east_border' = c(Inf)#, -140)
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
                                  DirPath=VAST_dir, 
                                  fine_scale = fine_scale,
                                  Save_Results=T )

# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, 
                      "knot_i"=Spatial_List$knot_i )

#plot data points and knots
plot(Spatial_List$latlon_i[,2:1])
points(Spatial_List$latlon_x[,2:1], col = 'red', pch = 16, cex = 1)

save.image(paste0(VAST_dir,'/Spatial_Settings_CrVa.RData') )
