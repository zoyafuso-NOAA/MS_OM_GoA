###############################################################################
## Project:      Univariate VAST model runs
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
##               Jim Thorson's VAST wiki example
##           (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:  Run single-species VAST models for 
##               Sebastes variabilis
##               Sebastes aleutus
##               Sebastes polyspinis
##               Run 5-fold Cross Validation for Each Model
###############################################################################

##################################################
####   Load packages 
##################################################
library(VAST)

##################################################
####   Set up model settings
####   Create new directories for each model
####   Set up directories
##################################################
rm(list = ls())

github_dir = 'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/'
VAST_dir = 'C:/Users/Zack Oyafuso/Google Drive/GOA_VAST_Runs/Single_Species/'

which_spp = c('Sebastes polyspinis', 
              'Sebastes variabilis',
              'Sebastes alutus')[3]

result_dir = paste0(VAST_dir, which_spp, '/')
if(!dir.exists(result_dir)) dir.create(result_dir)  

##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
master_data = read.csv(file = paste0(github_dir, 'data/GOA_multspp.csv') )

##################################################
####   Subset species
##################################################
data = subset(master_data, SPECIES_NAME == which_spp)
data$SPECIES_NAME = droplevels(data$SPECIES_NAME)

##################################################
####   Prepare the dataframe for catch-rate data in the VAST format
##################################################
Data_Geostat = data.frame( "spp" = data$SPECIES_NAME,
                           "Year" = data$YEAR,
                           "Catch_KG" = data$WEIGHT,
                           "AreaSwept_km2" = data$EFFORT,
                           "Vessel" = 0,
                           "Lat" = data$LATITUDE,
                           "Lon" = data$LONGITUDE)

##################################################
####   Assign 5 fold partitions of the data
##################################################
n_fold = 5
years = paste0(unique(Data_Geostat$Year))
NTime = length(unique(Data_Geostat$Year))

#Create unique stationID from the latlon. To make sure the ids are unique,
#we use the table function to make sure there are 7900 records (as of 2019)
Data_Geostat$latlon = paste0(Data_Geostat$Lat, Data_Geostat$Lon)
table(table(Data_Geostat$latlon))

#split Data_Geostat by year, then on each year-split, randomly assign 
#fold numbers to the each unique station
set.seed(2342)
foldno = lapply(X = split.data.frame(Data_Geostat, f = Data_Geostat$Year),
                FUN = function(test) {
                  unique_loc = unique(test$latlon)
                  fold_no = sample(x = 1:n_fold, 
                                   size = length(unique_loc), 
                                   replace = T)
                  return(split(unique_loc, fold_no))
                })

#Attach fold number to the Data_Geostat
for(iyear in years){
  for(ifold in paste(1:n_fold)){
    Data_Geostat[Data_Geostat$latlon %in% foldno[[iyear]][[ifold]] , 'fold'] =
      as.integer(ifold) 
  }
}

##################################################
####   Spatial settings: The following settings define the spatial resolution 
####   for the model, and whether to use a grid or mesh approximation
####   Stratification for results
##################################################
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
strata.limits <- data.frame('STRATA' = c("All_areas"),
                            'west_border' = c(-Inf),
                            'east_border' = c(Inf))

settings = FishStatsUtils::make_settings( 
  n_x = 350,   # Number of knots
  Region='Gulf_of_Alaska', 
  purpose="index2",
  strata.limits=strata.limits, 
  bias.correct=FALSE,
  FieldConfig =  c("Omega1"=1,  "Epsilon1"=1, "Omega2"=1,  "Epsilon2"=1),
  RhoConfig = c("Beta1"=0,   "Beta2"=0, "Epsilon1"=0,  "Epsilon2"=0),
  OverdispersionConfig = c("Eta1"=0,   "Eta2"=0), 
  ObsModel = c(2,0) ,
  Options = c("SD_site_density"=0, "SD_site_logdensity"=0, "Project_factors"=0),
  use_anisotropy = T)

##################################################
####   Fit the model and save output
##################################################
fit = fit_model( "settings"=settings,
                 "working_dir" = result_dir,
                 "Lat_i"=Data_Geostat[,'Lat'],
                 "Lon_i"=Data_Geostat[,'Lon'],
                 "t_i"=Data_Geostat[,'Year'],
                 "c_i"=as.numeric(Data_Geostat[,'spp'])-1,
                 "b_i"=Data_Geostat[,'Catch_KG'],
                 "a_i"=Data_Geostat[,'AreaSwept_km2'],
                 "v_i"=Data_Geostat[,'Vessel'],
                 "max_cells" = Inf,
                 "getJointPrecision"=TRUE,
                 "newtonsteps" = 1,
                 'test_fit' = F)
save(list = c('fit', 'Data_Geostat'), file = paste0(result_dir, '/fit.RData'))


##################################################
####   5-fold Cross Validation
##################################################
n_fold = 5
for(fI in 1:n_fold){ 
  if(!dir.exists(paste0(result_dir, 'CV_', fI))){
    dir.create(paste0(result_dir, 'CV_', fI))
  }
} 

# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Data_Geostat$fold == fI, TRUE, FALSE )
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model( "settings"=settings, 
                       "working_dir"=paste0(result_dir,'CV_', fI),
                       "Lat_i"=Data_Geostat[,'Lat'],
                       "Lon_i"=Data_Geostat[,'Lon'], 
                       "t_i"=Data_Geostat[,'Year'],
                       "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                       "b_i"=Data_Geostat[,'Catch_KG'],
                       "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                       "v_i"=Data_Geostat[,'Vessel'],
                       "PredTF_i"=PredTF_i, 
                       "Parameters"=fit$ParHat,
                       "getsd"=T,
                       "silent" = T,
                       "max_cells" = Inf,
                       "test_fit" = F,
                       "newtonsteps" = 1)
  
  # Save fit 
  save(list = 'fit_new',  file = paste0(result_dir,'CV_', fI, '/fit.RData'))
}

