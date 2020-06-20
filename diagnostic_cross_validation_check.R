######################################
## Diagnostics
######################################

library(VAST); library(RANN)

#Load Dataset
setwd("C:/Users/zack.oyafuso/Desktop/VAST_Runs/VAST_output6i/")

CV_df = data.frame(ifold = 1:5)
RRMSE = array(dim = c(5, ncol = 15, 11))

for(ifold in 1:5){
  load(paste0("CV_", ifold, '/fit.RData'))
  
  # #Final Gradient
  CV_df$max_grad[ifold] = max(abs(fit_new$parameter_estimates$diagnostics$final_gradient))
  
  #if hessian is positive definite, the sds are all be positive and the
  #eigenvalues of the covariance matrix are all positive
  sds = sqrt(diag(fit_new$parameter_estimates$SD$cov.fixed))
  CV_df$pdHess[ifold] = fit_new$parameter_estimates$SD$pdHess
  
  #check_fit chekcs bounds, TRUE is bad and FALSE is good
  CV_df$bound_check[ifold] = (check_fit(fit_new$parameter_estimates))
  
  withheld_idx = which(fit_new$data_list$PredTF_i == T)
  #Withheld data locations
  withheld_df = data.frame(
    idx =  withheld_idx,
    E_km = fit_new$spatial_list$loc_i[withheld_idx,'E_km'],
    N_km = fit_new$spatial_list$loc_i[withheld_idx,'N_km'],
    spp = 1+fit_new$data_frame$c_iz[withheld_idx],
    year = 1+fit_new$data_list$t_iz[withheld_idx,1],
    obs_density = (fit_new$data_frame$b_i/fit_new$data_frame$a_i)[withheld_idx]
  )
  
  
  #Grid locations
  loc_g = fit_new$spatial_list$loc_g
  
  #which grids are closest to each withheld data location 
  grid_idx = RANN::nn2(query = withheld_df[,c('E_km', 'N_km')], 
                       data = loc_g, k = 1)$nn.idx
  
  for(irow in 1:nrow(withheld_df)) {
    withheld_df$pred_density[irow] = fit_new$Report$D_gcy[grid_idx[irow],
                                                          withheld_df$spp[irow],
                                                          withheld_df$year[irow]] 
  }
  
  #RRMSE
  mean_pred_density = aggregate(pred_density ~ spp + year, data = withheld_df, 
                                FUN = mean)
  
  
  for(i in 1:15){
    for(itime in 1:11){
      split_df = subset(withheld_df, spp == i & year ==  unique(withheld_df$year)[itime])
      temp_RMSE = sqrt(mean((split_df$obs_density - split_df$pred_density)^2))
      temp_mean_pred_density = mean_pred_density[mean_pred_density$spp == i &
                                                   mean_pred_density$year == unique(withheld_df$year)[itime],
                                                 'pred_density']
      RRMSE[ifold, i, itime] = temp_RMSE / temp_mean_pred_density
    }
    
  }
  print(paste0('Done with fold ', ifold))
}

CV_df
round(apply(RRMSE, MARGIN = 1:2, mean), 2)
