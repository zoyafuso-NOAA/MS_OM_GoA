###############################################################################
## Project:      VAST model runs with Cross-Validation
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
##               Jim Thorson's VAST wiki example
##           (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:  Run multispecies VAST model with different settings wrt
##               number of factors, species set, density covariates, etc.
###############################################################################

##################################################
####   Load packages 
##################################################
library(VAST)
library(TMBhelper)

##################################################
####   Set up model settings
####   Create new directories for each model
####   Set up directories
##################################################
rm(list = ls())
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2)[1]

model_settings = data.frame(factorno = 2,
                            modelno = paste0(9, letters[1:2]),
                            stringsAsFactors = F)

irow = 1 #Specify which model 
factorno = model_settings$factorno[irow]
modelno = model_settings$modelno[irow]

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', modelno, '/')

if(!dir.exists(VAST_dir)) dir.create(VAST_dir)

##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
data = read.csv(file = paste0(github_dir, 'data/GOA_multspp.csv') )
spp_df = read.csv(paste0(github_dir, "data/spp_df.csv"), 
                  check.names=F, header = T, row.names = 'modelno')

##################################################
####   Prepare the dataframe for catch-rate data in the VAST format
##################################################
Data_Geostat = data.frame( "spp" = data$SPECIES_NAME,
                           "Year" = data$YEAR,
                           "Catch_KG" = data$WEIGHT,
                           "AreaSwept_km2" = data$EFFORT,
                           "Vessel" = 0,
                           "Lat" = data$LATITUDE,
                           "Lon" = data$LONGITUDE,
                           #Centered depth variables
                           "LOG_DEPTH" = data$DEPTH,
                           "LOG_DEPTH2" = data$DEPTH^2 )

##################################################
#### Subset Species and drop factor levels of unused Species
##################################################
which_spp = unlist(spp_df[modelno,])

Data_Geostat = subset(Data_Geostat, spp %in% names(which_spp)[which_spp])
Data_Geostat$spp = droplevels(Data_Geostat$spp)

##################################################
####   Assign 5 fold partitions of the data
##################################################
n_fold = 5
years = paste0(unique(Data_Geostat$Year))
NTime = length(unique(Data_Geostat$Year))

#Sort Data_Geostat
Data_Geostat = Data_Geostat[order(Data_Geostat$Year, Data_Geostat$spp),]

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
####   Save Spatial Object
##################################################
save(list = c('Data_Geostat'),
     file = paste0(VAST_dir,'Spatial_Settings_CrVa.RData') )
# load(paste0(VAST_dir, '/Spatial_Settings_CrVa.RData'))

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
  FieldConfig =  c("Omega1"=factorno,  "Epsilon1"=factorno, 
                   "Omega2"=factorno,  "Epsilon2"=factorno),
  RhoConfig = c("Beta1"=0,   "Beta2"=0, 
                "Epsilon1"=0,  "Epsilon2"=0),
  OverdispersionConfig = c("Eta1"=0,   "Eta2"=0), 
  ObsModel = c(2,0) ,
  Options = c("SD_site_density"=0, "SD_site_logdensity"=0, "Project_factors"=0),
  use_anisotropy = T)

##################################################
####   Import "true" and not interpolated covariate 
####   data if using depth covariates
##################################################
load( paste0(github_dir, 'data/Extrapolation_depths.RData'))

n_g = nrow(Extrapolation_depths) #number of grid cells
n_t = diff(range(Data_Geostat$Year)) + 1 #Number of total years
n_p = 2 #two density covariates

X_gtp = array(dim = c(n_g, n_t, n_p) )
for(i in 1:Save$TmbData$n_t) {
  X_gtp[,i,] = as.matrix(Extrapolation_depths[,c('DEPTH', 'DEPTH2')])
}

##################################################
####   Fit the model and save output
####   If adding density covariates, the VAST version that I am using for MRAN
####   will give an error if you supply it with X_gtp. So fit_model_X_GTP.R
####   modifies the X_gtp object after make_data() is called in fit_model. 
##################################################
source(paste0(github_dir, 'fit_model_X_GTP.R'))
fit = fit_model( "settings"=settings,
                 "working_dir" = VAST_dir,
                 "Lat_i"=Data_Geostat[,'Lat'],
                 "Lon_i"=Data_Geostat[,'Lon'],
                 "t_i"=Data_Geostat[,'Year'],
                 "c_i"=as.numeric(Data_Geostat[,'spp'])-1,
                 "b_i"=Data_Geostat[,'Catch_KG'],
                 "a_i"=Data_Geostat[,'AreaSwept_km2'],
                 "v_i"=Data_Geostat[,'Vessel'],
                 # "formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                 # "covariate_data" = cbind(Data_Geostat[,c('Lat', 'Lon',
                 #                                          'LOG_DEPTH',
                 #                                          'LOG_DEPTH2',
                 #                                          'Catch_KG')],
                 # Year = NA),
                 "max_cells" = Inf,
                 "getJointPrecision"=TRUE,
                 "newtonsteps" = 1#,
                 # 'X_gtp' = X_gtp
)
ParHat = fit$ParHat

## Save Fit
save(list = 'fit', file = paste0(VAST_dir, '/fit.RData'))

# load(paste0(VAST_dir, 'VAST_MS_GoA_Run.RData'))
# ParHat = Save$ParHat

# load(paste0(VAST_dir, 'fit.RData'))
# ParHat = fit$ParHat

##################################################
####   5-fold Cross Validation
##################################################
n_fold = 5
for(fI in 1:n_fold){ 
  if(!dir.exists(paste0(VAST_dir, 'CV_', fI))){
    dir.create(paste0(VAST_dir, 'CV_', fI))
  }
} 

# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Data_Geostat$fold == fI, TRUE, FALSE )
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model( "settings"=settings, 
                       "working_dir"=paste0(VAST_dir,'CV_', fI),
                       "Lat_i"=Data_Geostat[,'Lat'],
                       "Lon_i"=Data_Geostat[,'Lon'], 
                       "t_i"=Data_Geostat[,'Year'],
                       "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                       "b_i"=Data_Geostat[,'Catch_KG'],
                       "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                       "v_i"=Data_Geostat[,'Vessel'],
                       "PredTF_i"=PredTF_i, 
                       "Parameters"=ParHat,
                       "getsd"=T,
                       "silent" = T,
                       "formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                       "covariate_data" = cbind(Data_Geostat[,c('Lat', 'Lon',
                                                                'LOG_DEPTH',
                                                                'LOG_DEPTH2',
                                                                'Catch_KG')],
                                                Year = NA),
                       "max_cells" = Inf,
                       "newtonsteps" = 1)#,
  # 'X_gtp' = X_gtp)
  
  # Save fit 
  save(list = 'fit_new',  file = paste0(VAST_dir,'CV_', fI, '/fit.RData'))
}

