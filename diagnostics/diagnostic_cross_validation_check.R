###############################################################################
## Project:       Relative Root Mean Square Error Calculations
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate RRMSE for each 5-fold CrossValidation Run
###############################################################################
rm(list = ls())

##################################################
####   Import Required packages
##################################################
library(VAST)
library(RANN)

##################################################
####   Set up directories
##################################################
VAST_dir = "C:/Users/zack.oyafuso/Desktop/VAST_Runs/VAST_output6f/"

##################################################
####   Result Objects
##################################################
CV_df = data.frame(ifold = 1:5)
RRMSE = array(dim = c(5, ncol = 15, 11))

for (ifold in 1:5){
  
  ##################################################
  ####   Load Result Object
  ##################################################
  load(paste0(VAST_dir, "CV_", ifold, '/fit.RData'))
  
  ##################################################
  ####   Some Convergence/Hessian Checks
  ##################################################
  #Final maximum absolute gradient
  CV_df$max_grad[ifold] <-
    max(abs(fit_new$parameter_estimates$diagnostics$final_gradient))
  
  #Is hessian is positive definite?
  CV_df$pdHess[ifold] <- fit_new$parameter_estimates$SD$pdHess
  
  #check_fit chekcs bounds, TRUE is bad and FALSE is good
  CV_df$bound_check[ifold] <- check_fit(fit_new$parameter_estimates)
  
  ##################################################
  ####   Calculate RRMSE
  ##################################################
  withheld_idx <- which(fit_new$data_list$PredTF_i == T)
  
  #Withheld data locations
  withheld_df <- data.frame(
    idx =  withheld_idx,
    E_km = fit_new$spatial_list$loc_i[withheld_idx,'E_km'],
    N_km = fit_new$spatial_list$loc_i[withheld_idx,'N_km'],
    spp = 1+fit_new$data_frame$c_iz[withheld_idx],
    year = 1+fit_new$data_list$t_iz[withheld_idx,1],
    obs_density = (fit_new$data_frame$b_i/fit_new$data_frame$a_i)[withheld_idx]
  )
  
  #Grid locations
  loc_g <- fit_new$spatial_list$loc_g
  
  #Year indices
  years = unique(withheld_df$year)
  
  #which grids are closest to each withheld data location
  grid_idx <- RANN::nn2(query = withheld_df[,c('E_km', 'N_km')],
                        data = loc_g,
                        k = 1)$nn.idx
  
  #For each station, locate cell/spp/year density predictions that correspond
  #to the location of the data
  for (irow in 1:nrow(withheld_df)) {
    withheld_df$pred_density[irow] <-
      fit_new$Report$D_gcy[grid_idx[irow],
                           withheld_df$spp[irow],
                           withheld_df$year[irow]]
  }
  
  #Calculate the mean observed density for the calculation of the RRMSE
  #There are many ways to do this and I chose for the normalizing constant to
  #be specific to year and species.
  mean_pred_density = aggregate(obs_density ~ spp + year, data = withheld_df,
                                FUN = mean)
  
  #Calculate RRMSE for each year for each species using the mean observed
  #density calculated above
  for (ispp in 1:15) {
    for (itime in 1:11) {
      split_df <- subset(withheld_df,
                         spp == ispp & year == years[itime])
      temp_RMSE <- sqrt(mean((split_df$obs_density - split_df$pred_density)^2))
      temp_mean_pred_density <-
        mean_pred_density[mean_pred_density$spp == ispp
                          & mean_pred_density$year == years[itime],
                          'obs_density']
      RRMSE[ifold, ispp, itime] = temp_RMSE / temp_mean_pred_density
    }
    
  }
  print(paste0('Done with fold ', ifold))
}

CV_df

#Average RRMSE across years
round(apply(RRMSE, MARGIN = 1:2, mean), 2)
