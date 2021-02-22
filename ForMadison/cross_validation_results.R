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
library(tidyr)

##################################################
####   Set up directores
##################################################
VAST_dir <-  "G:/Oyafuso/VAST_Runs_EFH/Single_Species/" 
github_dir <- "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/"

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, "data/Extrapolation_depths.RData"))

which_spp <- sort(c("Sebastes polyspinis", "Sebastes variabilis", 
                    "Sebastes brevispinis", "Microstomus pacificus",
                    "Lepidopsetta polyxystra", "Lepidopsetta bilineata",
                    "Hippoglossus stenolepis", "Hippoglossoides elassodon",
                    "Glyptocephalus zachirus", "Gadus macrocephalus",
                    "Gadus chalcogrammus", "Sebastes B_R", "Sebastes alutus",
                    "Atheresthes stomias", "Sebastolobus alascanus",
                    "Anoplopoma fimbria", "Beringraja spp", "Octopus spp",
                    "Pleurogrammus monopterygius", "Sebastes borealis",
                    "Sebastes variegatus", "Squalus suckleyi"))

ns <- length(which_spp)
n_folds <- 10
n_years <- 11

##################################################
####   Result Objects
##################################################
cv_df <- expand.grid(species = which_spp, 
                     depth = c(T, F),
                     fold = 1:n_folds,
                     stringsAsFactors = F)

cv_df[,c("max_grad", "pdHess", "bound_check", "pred_nll", 
         "RMSE", "RRMSE", "MAE", "RMAE")] <- NA

for (irow in (1:nrow(cv_df)) ) {
  
  #Load fitted object
  result_dir <- paste0(VAST_dir, cv_df$species[irow], 
                       ifelse(cv_df$depth[irow],  "_depth", ""), "/")
  filename <- paste0(result_dir, "CV_", cv_df$fold[irow], "/fit.RData")
  
  if ( file.exists(filename) ){
    
    ## Load crossvalidation result and assign temporary result objects
    load(filename)
    pars <- fit_new$parameter_estimates
    spatial_list <- fit_new$spatial_list
    temp_df <- fit_new$data_frame
    
    #Final Gradient
    cv_df$max_grad[irow] <- max(abs(pars$diagnostics$final_gradient))
    
    #Check whether hessian matrix is positive definite
    cv_df$pdHess[irow] <- pars$SD$pdHess
    
    #check_fit chekcs bounds, TRUE is bad and FALSE is good
    cv_df$bound_check[irow] <- check_fit(pars)
    
    cv_df$pred_nll[irow] <- fit_new$Report$pred_jnll
    
    #Extract incides of withheld data
    withheld_idx <- which(fit_new$data_list$PredTF_i == T)
    
    #Withheld data locations, species/year indices, observed CPUE
    withheld_df <- data.frame(
      idx =  withheld_idx,
      E_km = spatial_list$loc_i[withheld_idx, "E_km"],
      N_km = spatial_list$loc_i[withheld_idx, "N_km"],
      year = 1 + fit_new$data_list$t_i[withheld_idx],
      obs_density = (temp_df$b_i / temp_df$a_i)[withheld_idx]
    )
    
    #Extrapolation Grid locations
    loc_g <- spatial_list$loc_g
    
    #which extrapoaltion grid cells are closest to each withheld data location 
    grid_idx <- RANN::nn2(query = withheld_df[,c("E_km", "N_km")], 
                          data = loc_g, 
                          k = 1)$nn.idx
    
    # for each withheld datum, assign its predicted density to that of its
    # nearest neighboring extrapolation grid cell
    for (jrow in 1:nrow(withheld_df)) {
      withheld_df$pred_density[jrow] <- 
        fit_new$Report$D_gct[grid_idx[jrow], , withheld_df$year[jrow]]
    }
    
    #Calculate mean absolute error and root mean square error
    cv_df$MAE[irow] = mean(abs(withheld_df$obs_density - 
                                 withheld_df$pred_density))
    
    cv_df$RMSE[irow] = sqrt(mean((withheld_df$obs_density - 
                                    withheld_df$pred_density)^2))
    
    #Calculate mean obs density for calculation of RRMSE and RMAE
    mean_obs_density <- mean(temp_df$b_i / temp_df$a_i)
    cv_df$RRMSE[irow] <- cv_df$RMSE[irow] / mean_obs_density
    cv_df$RMAE[irow] <- cv_df$MAE[irow] / mean_obs_density
    
    #Update progress
    print(paste0("Done with: ", cv_df$species[irow], ", ", 
                 ifelse(cv_df$depth[irow], "Depth, ", "No Depth, "),
                 "Fold Number ", cv_df$fold[irow]))
  }
}

##################################################
####  Calculate summed predicted NLL across folds and Mean RRMSE across folds
##################################################
tidyr::spread(data = aggregate(pred_nll ~ species + depth,
                               data = cv_df,
                               FUN = sum),
              key = "depth",
              value = "pred_nll")

RMAE <- tidyr::spread(data = aggregate(RMAE ~ species + depth,
                                       data = cv_df,
                                       FUN = mean),
                      key = "depth",
                      value = "RMAE")

RMAE$depth_in_model <- c(F, T)[apply(X = RMAE[,-1],
                                     MARGIN = 1,
                                     FUN = which.min)]

RRMSE <- tidyr::spread(data = aggregate(RRMSE ~ species + depth,
                                        data = cv_df,
                                        FUN = mean),
                       key = "depth",
                       value = "RRMSE")

RRMSE$depth_in_model <- c(F, T)[apply(X = RRMSE[,-1],
                                      MARGIN = 1,
                                      FUN = which.min)]

prednll <- tidyr::spread(data = aggregate(pred_nll ~ species + depth,
                                       data = cv_df,
                                       FUN = mean),
                      key = "depth",
                      value = "pred_nll")

prednll$depth_in_model <- c(F, T)[apply(X = prednll[,-1],
                                     MARGIN = 1,
                                     FUN = which.min)]

##################################################
####   Create the result object that would go into the optimizations
##################################################
N <- nrow(Extrapolation_depths)
D_gct = Index <- array(dim = c(N, ns, 24), 
                       dimnames = list(NULL, which_spp, NULL))

for(ispp in 1:ns){
  depth_in_model <- RRMSE$depth_in_model[ispp]
  
  result_dir = paste0(VAST_dir, RMSE$species[ispp],
                      ifelse(depth_in_model, "_depth", ""),
                      "/")
  
  load(paste0(result_dir, "/fit.RData"))
  
  D_gct[,ispp,] = fit$Report$D_gct[,1,]
  Index[,ispp,] = fit$Report$Index_gctl[,1, ,1]
  
  print(result_dir)
}

##################################################
####   Save
##################################################
save("RRMSE", file = paste0(github_dir, "/data/RRMSE_VAST_models.RData"))
save("D_gct", file = paste0(github_dir, "/data/fit_density.RData") )
save("Index", file = paste0(github_dir, "/data/fit_index.RData") )

