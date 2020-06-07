######################################
## Diagnostics
######################################

library(VAST); library(RANN)

#Load Dataset
setwd("C:/Users/zack.oyafuso/Desktop/VAST_Runs/VAST_output6i/")

CV_df = data.frame(ifold = 1:5)
RRMSE = matrix(nrow = 5, ncol = 15)

for(ifold in 1:4){
  load(paste0("CV_", ifold, '/fit.RData'))
  
  # #Final Gradient
  CV_df$max_grad[ifold] = max(abs(fit_new$parameter_estimates$diagnostics$final_gradient))
  
  #if hessian is positive definite, the sds are all be positive and the
  #eigenvalues of the covariance matrix are all positive
  sds = sqrt(diag(fit_new$parameter_estimates$SD$cov.fixed))
  CV_df$all_sd_positive[ifold] = all(sds > 0)
  CV_df$all_eigen_positive[ifold] = all(eigen(fit_new$parameter_estimates$SD$cov.fixed)$values>0)
  
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
  grid_idx = RANN::nn2(query = withheld_df[,c('E_km', 'N_km')], data = loc_g, k = 1)$nn.idx
  
  for(irow in 1:nrow(withheld_df)) {
    withheld_df$pred_density[irow] = fit_new$Report$D_gcy[grid_idx[irow],
                                                          withheld_df$spp[irow],
                                                          withheld_df$year[irow]] 
  }
  
  #RRMSE
  mean_pred_density = aggregate(pred_density ~ spp, data = withheld_df, 
                                FUN = mean)$pred_density
  split_df = split.data.frame(withheld_df, f = withheld_df$spp)
  
  for(i in 1:15){
    RRMSE[ifold, i] = with(split_df[[i]], 
                             sqrt(mean((obs_density-pred_density)^2))/mean_pred_density[i] )
  }
  
  
  
}

colMeans(RRMSE, na.rm = T)
