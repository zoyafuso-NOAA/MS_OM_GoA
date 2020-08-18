###############################################################################
## Project:       Cross-Validation Results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For each CV run, calculate relative root mean squre error of
##                density predictions
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(VAST)
library(RANN)

##################################################
####   Set up directores
##################################################
VAST_dir = 'C:/Users/Zack Oyafuso/Google Drive/GOA_VAST_Runs/Single_Species/'

which_spp = c('Sebastes polyspinis', 
              'Sebastes variabilis',
              'Sebastes alutus')[2]

result_dir = paste0(VAST_dir, which_spp, '/')

##################################################
####   Result Objects
##################################################
nFold = 5
nTime = 11
CV_df = data.frame(ifold = 1:nFold)
RRMSE = array(dim = c(nFold, nTime))

for(ifold in 1:nFold){
  load(paste0(result_dir, "CV_", ifold, '/fit.RData'))
  
  #Final Gradient
  CV_df$max_grad[ifold] = max(abs(
    fit_new$parameter_estimates$diagnostics$final_gradient))
  
  #Check whether hessian matrix is positive definite
  CV_df$pdHess[ifold] = fit_new$parameter_estimates$SD$pdHess
  
  #check_fit chekcs bounds, TRUE is bad and FALSE is good
  CV_df$bound_check[ifold] = (check_fit(fit_new$parameter_estimates))
  
  #Extract incides of withheld data
  withheld_idx = which(fit_new$data_list$PredTF_i == T)
  
  #Withheld data locations, species/year indices, observed CPUE
  withheld_df = data.frame(
    idx =  withheld_idx,
    E_km = fit_new$spatial_list$loc_i[withheld_idx,'E_km'],
    N_km = fit_new$spatial_list$loc_i[withheld_idx,'N_km'],
    spp = 1+fit_new$data_frame$c_iz[withheld_idx],
    year = 1+fit_new$data_list$t_iz[withheld_idx,1],
    obs_density = (fit_new$data_frame$b_i/fit_new$data_frame$a_i)[withheld_idx]
  )
  
  #Extrapolation Grid locations
  loc_g = fit_new$spatial_list$loc_g
  
  #which extrapoaltion grid cells are closest to each withheld data location 
  grid_idx = RANN::nn2(query = withheld_df[,c('E_km', 'N_km')], 
                       data = loc_g, 
                       k = 1)$nn.idx
  
  for(irow in 1:nrow(withheld_df)) {
    temp_density = fit_new$Report$D_gcy[grid_idx[irow],
                                        withheld_df$spp[irow],
                                        withheld_df$year[irow]]
    withheld_df$pred_density[irow] =  temp_density
  }
  
  #Calculate mean predicted density for calculation of RRMSE
  mean_pred_density = aggregate(obs_density ~ year, data = withheld_df, 
                                FUN = mean)
  
  for(itime in 1:nTime){
    split_df = subset(withheld_df, year ==  unique(withheld_df$year)[itime])
    temp_RMSE = sqrt(mean((split_df$obs_density - split_df$pred_density)^2))
    temp_mean_pred_density = mean_pred_density[itime, 'obs_density']
    RRMSE[ifold, itime] = temp_RMSE / temp_mean_pred_density
  }
  
  print(paste0('Done with fold ', ifold))
}

CV_df
rowMeans(RRMSE)
mean(RRMSE)
