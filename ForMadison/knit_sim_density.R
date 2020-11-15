###############################################################################
## Project:       
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/Single_Species/"

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, '/data/optimization_data.RData'))
load(paste0(github_dir, '/data/RMSE_VAST_models.RData'))
load(paste0(github_dir, '/data/fit_density.RData'))

##################################################
####   Knit together density predictions from the original fits and the
####   the ten sets of density predictions
##################################################
pred_density <- list(pred_density = D_gct[,,Years2Include],
                     plus_obs_error = array(dim = c(10, N, ns, NTime)) )

for (ispp in 1:ns) {
  chosen_model = paste0(VAST_dir,
                        sci_names[ispp],
                        ifelse(RMSE$depth_in_model[ispp],
                               "_depth",
                               ""), "/")
  for (iter in 1:10) {
    #Load fitted model
    load(paste0(chosen_model, "sim_runs/iter_", iter, "/fit.RData"))
    pred_density$plus_obs_error[iter, ,ispp ,] <- fit$Report$D_gct[, , Years2Include]
    
    print(paste("Finished with", sci_names[ispp], "and iteration", iter))
  }
}

##################################################
#Save Predicted Density, move up
##################################################
save(list = "pred_density",
     file = paste0(dirname(VAST_dir), "/sim_density.RData"))
