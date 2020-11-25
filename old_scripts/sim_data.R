###############################################################################
## Project:       Simulate Data 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   simulate 10 sets of data for each species and save 
###############################################################################
rm(list = ls())

##################################################
####    Import required packages
##################################################
library(FishStatsUtils)

##################################################
####   Set up directories
####
####   Set up some constants of the optimization
####   Multispeceis: Spatiotemporal Variance, species specific CV constraints
####   Single_Species: Spatiotemporal Variance, univariate optimization, 
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine],
                     "GitHub/Optimal_Allocation_GoA/")

VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/"

##################################################
####    Load Data
##################################################
load(paste0(github_dir, "data/RMSE_VAST_models.RData"))

##################################################
####    Result Object
##################################################
sim_data <- array(dim = c(10, 7900, 15))


for (ispp in 1:ns) {
  #Load VAST dll
  dll_name <- paste0(VAST_dir, 
                     "Single_Species/",
                     sci_names[ispp],
                     ifelse(RMSE[ispp, "depth_in_model"], "_depth", ""), "/" )

  dyn.load(paste0(dll_name, "/VAST_v12_0_0.dll"))

  #Load fitted model
  load( paste0(dll_name, "fit.RData") )
  
  for (iter in 1:10) {
    temp_sim <- FishStatsUtils::simulate_data(fit = fit,
                                              type = 1,
                                              random_seed = iter*1000)
    
    sim_data[iter, , ispp] = temp_sim$b_i
    
    print(paste("Finished with", sci_names[ispp], "Iteration", iter ))
  }
  
  dyn.unload(paste0(dll_name, "/VAST_v12_0_0.dll"))
}

##################################################
####    Save Object
##################################################
save(list = "sim_data", 
     file = paste0(VAST_dir, "sim_data.RData"))
