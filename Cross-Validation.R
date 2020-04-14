########################################
## Cross Validation
## Number of Species Factors, 2-5
########################################

# Load packages
library(TMB)
library(VAST)

###################################
## Set up directories
###################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2, 'VM' = 3)[2]

setwd(paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output6j', '/'))

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/')

VAST_modelno = data.frame(factorn = c(2,3,5),
                     modelno = paste0(6, c('f', 'd', 'h')),
                     stringsAsFactors = F)
irow=1

VAST_model = VAST_modelno$modelno[irow]
factorno = VAST_modelno$factorn[irow]

if(!dir.exists(paste0('Factor_', factorno)))
  dir.create(paste0('Factor_', factorno))

###############################
## Load Model Information, Fitted Parameters
###############################
#Load Model SEttings
load(paste0(VAST_dir, 'VAST_output', VAST_model, "/Model_Settings.RData"))

load(paste0(VAST_dir, 'VAST_output6j/Spatial_Settings_CrVa.RData'))
settings = make_settings( n_x=n_x, 
                          Region='Gulf_of_Alaska', 
                          purpose="index",
                          strata.limits=strata.limits, 
                          bias.correct=FALSE,
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig, 
                          ObsModel = ObsModel,
                          Options = Options,
                          use_anisotropy = T)

# Fit the model and a first time and record MLE
fit = fit_model( "settings"=settings,
                 "Lat_i"=Data_Geostat[,'Lat'],
                 "Lon_i"=Data_Geostat[,'Lon'],
                 "t_i"=Data_Geostat[,'Year'],
                 "c_i"=as.numeric(Data_Geostat[,'spp'])-1,
                 "b_i"=Data_Geostat[,'Catch_KG'],
                 "a_i"=Data_Geostat[,'AreaSwept_km2'],
                 "v_i"=Data_Geostat[,'Vessel'] )
ParHat = fit$ParHat
save(list = 'fit', file = paste0(VAST_dir, '/Factor_', factorno, '/fit.RData'))

# Generate partitions in data
prednll_f = rep(NA, n_fold )

# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Data_Geostat$fold == fI, TRUE, FALSE )
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model( "settings"=settings, 
                       "Lat_i"=Data_Geostat[,'Lat'],
                       "Lon_i"=Data_Geostat[,'Lon'], 
                       "t_i"=Data_Geostat[,'Year'],
                       "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                       "b_i"=Data_Geostat[,'Catch_KG'],
                       "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                       "v_i"=Data_Geostat[,'Vessel'],
                       "PredTF_i"=PredTF_i, 
                       "Parameters"=ParHat, 
                       "getsd"=FALSE )
  
  
  #Finally, we bundle and save output
  
  # Save fit to out-of-bag data
 prednll_f[fI] = fit_new$Report$pred_jnll
 save(list = 'prednll_f', 
      file = paste0(VAST_dir, '/Factor_', factorno, '/CrVa_Run.RData'))
}

# Check fit to all out=of-bag data and use as metric of out-of-bag performance
sum( prednll_f )