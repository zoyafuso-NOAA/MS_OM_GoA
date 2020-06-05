######################################
## Diagnostics
######################################

library(VAST); library(RANN)

#Load Dataset
setwd("C:/Users/zack.oyafuso/Desktop/VAST_Runs/VAST_output6f/")

CV_df = data.frame(ifold = 1:5)

for(ifold in 1:5){
  load(paste0("CV_", ifold, '/fit.RData'))
  
  # #Final Gradient
  CV_df$max_grad[ifold] = max(abs(fit_new$parameter_estimates$diagnostics$final_gradient))
  
  
  #if hessian is positive definite, the sds are all be positive and the
  #eigenvalues of the covariance matrix are all positive
  CV_df$all_sd_positive[ifold] = all(df$sd > 0)
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
  grid_idx = RANN::nn2(query = withheld_loc, data = loc_g, k = 1)$nn.idx
  
  for(irow in 1:nrow(withheld_df)) {
    withheld_df$pred_density[irow] = fit_new$Report$D_gcy[grid_idx[irow],
                                                          withheld_df$spp[irow],
                                                          withheld_df$year[irow]] 
  }
  
  CV_df$RRMSE[ifold] = with(withheld_df, 
                            sqrt(mean((obs_density-pred_density)^2))/mean(pred_density) )
  
}
