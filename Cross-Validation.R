########################################
## Cross Validation
## Number of Species Factors, 2-5
########################################
rm(list = ls())

# Load packages
library(VAST)
library(TMBhelper)
# library(devtools)
# install_local('C:/Users/zack.oyafuso/Downloads/FishStatsUtils-2.6.0/',
#               force = T)

###################################
## Set up directories
###################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2, 'VM' = 3)[1]

setwd(paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
               'C:/Users/zack.oyafuso/Desktop/',
               'C:/Users/zack.oyafuso/Desktop/')[which_machine],
             'VAST_Runs/VAST_output7', '/'))

model_settings = data.frame(factorno = 2:6,
                            modelno = paste0(7, letters[1:5]),
                            stringsAsFactors = F)

irow = 1
factorno = model_settings$factorno[irow]
modelno = model_settings$modelno[irow]

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output7/')


if(!dir.exists(paste0('VAST_output', modelno, '/'))){
  dir.create(paste0('VAST_output', modelno, '/'))
}


###############################
## Spatial settings: The following settings define the spatial resolution 
## for the model, and whether to use a grid or mesh approximation
## Stratification for results
###############################
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
n_x = 350   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame(
  'STRATA' = c("All_areas"),#, "west_of_140W"),
  'west_border' = c(-Inf),#, -Inf),
  'east_border' = c(Inf)#, -140)
)

load(paste0(VAST_dir, 'Spatial_Settings_CrVa.RData'))

settings = make_settings( n_x=n_x, 
                          Region='Gulf_of_Alaska', 
                          purpose="index2",
                          strata.limits=strata.limits, 
                          bias.correct=FALSE,
                          FieldConfig =  c("Omega1"=factorno, 
                                           "Epsilon1"=factorno, 
                                           "Omega2"=factorno, 
                                           "Epsilon2"=factorno) ,
                          RhoConfig = c("Beta1"=0, "Beta2"=0, 
                                        "Epsilon1"=0, "Epsilon2"=0),
                          OverdispersionConfig = c("Eta1"=0, "Eta2"=0), 
                          ObsModel = c(2,0) ,
                          Options = c("SD_site_density"=0, 
                                      "SD_site_logdensity"=0, 
                                      "Project_factors" = 0),
                          use_anisotropy = T)

# Fit the model and a first time and record MLE
fit = fit_model( "settings"=settings,
                 "working_dir" = paste0(VAST_dir, 'VAST_output', modelno, '/'),
                 "Lat_i"=Data_Geostat[,'Lat'],
                 "Lon_i"=Data_Geostat[,'Lon'],
                 "t_i"=Data_Geostat[,'Year'],
                 "c_i"=as.numeric(Data_Geostat[,'spp'])-1,
                 "b_i"=Data_Geostat[,'Catch_KG'],
                 "a_i"=Data_Geostat[,'AreaSwept_km2'],
                 "v_i"=Data_Geostat[,'Vessel'],
                 "formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                 "covariate_data" = cbind(Data_Geostat[,c('Lat', 'Lon', 
                                                          'LOG_DEPTH',
                                                          'LOG_DEPTH2',
                                                          'Catch_KG')], 
                                          Year = NA),
                 "max_cells" = Inf,
                 "getJointPrecision"=TRUE,
                 "newtonsteps" = 1
)
ParHat = fit$ParHat

###############################
## Save Fit
###############################
save(list = 'fit', file = paste0(VAST_dir, 'VAST_output', modelno, '/fit.RData'))


###################################
## 10-fold Cross Validation
###################################
n_fold = 10
for(fI in 1:n_fold){ 
  if(!dir.exists(paste0(VAST_dir,'VAST_output', modelno,'/CV_', fI))){
    dir.create(paste0(VAST_dir,'VAST_output', modelno,'/CV_', fI))}
  } 

# Generate partitions in data
prednll_f = rep(NA, n_fold )

# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Data_Geostat$fold == fI, TRUE, FALSE )
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model( "settings"=settings, 
                       "working_dir"=paste0(VAST_dir,'VAST_output',
                                            modelno,'/CV_', fI),
                       "Lat_i"=Data_Geostat[,'Lat'],
                       "Lon_i"=Data_Geostat[,'Lon'], 
                       "t_i"=Data_Geostat[,'Year'],
                       "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                       "b_i"=Data_Geostat[,'Catch_KG'],
                       "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                       "v_i"=Data_Geostat[,'Vessel'],
                       "PredTF_i"=PredTF_i, 
                       "Parameters"=ParHat,
                       "getsd"=F,
                       "formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                       "covariate_data" = cbind(Data_Geostat[,c('Lat', 'Lon', 
                                                                'LOG_DEPTH',
                                                                'LOG_DEPTH2',
                                                                'Catch_KG')], 
                                                Year = NA),
                       "max_cells" = Inf,
                       "newtonsteps" = 1)
  
  
  #Finally, we bundle and save output
  
  # Save fit to out-of-bag data
  prednll_f[fI] = fit_new$Report$pred_jnll
  save(list = 'prednll_f', 
       file = paste0(VAST_dir, 'VAST_output', modelno, '/CrVa_Run.RData'))
}

# Check fit to all out=of-bag data and use as metric of out-of-bag performance
sum( prednll_f )

