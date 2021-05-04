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
# VAST_dir <-  "G:/Oyafuso/VAST_Runs_EFH/Single_Species/" 
# VAST_data_dir <-  "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/" 
# github_dir <- "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/"

VAST_dir <-  "C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Single_Species/" 
VAST_data_dir <-  "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/" 
github_dir <- "C:/Users/Zack Oyafuso/Documents/GitHub/Optimal_Allocation_GoA/"
##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, "data/Extrapolation_depths.RData"))
goa_data <- read.csv(paste0(VAST_data_dir, "data/GOA_multspp.csv"))
which_spp <- gsub(x = grep(x = dir(VAST_dir), pattern = "_depth", value = TRUE),
                  pattern = "_depth",
                  replacement = "")

ns <- length(which_spp)
n_folds <- 10
n_years <- length(unique(goa_data$YEAR))
which_years <- min(unique(goa_data$YEAR)):max(unique(goa_data$YEAR)) 
which_years <- which_years %in% unique(goa_data$YEAR)

##################################################
####   Result Objects
##################################################
cv_df <- expand.grid(species = which_spp, 
                     depth = c(T, F),
                     fold = 1:n_folds,
                     stringsAsFactors = F)

cv_df[,c("max_grad", "pdHess", "pred_nll", "RMSE")] <- NA

# for (irow in seq(length.out = nrow(cv_df))[is.na(cv_df$pdHess)] )  {
for (irow in 1:nrow(cv_df) )  {
  # for (irow in 132)  {
  #Load fitted object
  result_dir <- paste0(VAST_dir, cv_df$species[irow], 
                       ifelse(cv_df$depth[irow],  "_depth", ""), "/")
  filename <- paste0(result_dir, "CV_", cv_df$fold[irow], "/fit.RData")
  
  if ( file.exists(filename) ){
    
    ## Load crossvalidation result and assign temporary result objects
    load(filename)
    pars <- fit_new$parameter_estimates
    temp_df <- fit_new$data_frame
    
    #Final Gradient
    cv_df$max_grad[irow] <- max(abs(pars$diagnostics$final_gradient))
    
    #Check whether hessian matrix is positive definite
    cv_df$pdHess[irow] <- pars$SD$pdHess
    
    #Extract the out-of-bag predicted NLL
    cv_df$pred_nll[irow] <- fit_new$Report$pred_jnll
    
    #Extract observed and predicted cpue of the withheld data
    withheld_idx <- which(fit_new$data_list$PredTF_i == T)
    
    obs_cpue <- (temp_df$b_i / temp_df$a_i)[withheld_idx]
    pred_cpue <- fit_new$Report$D_i[withheld_idx]
    
    #Calculate mean absolute error and root mean square error
    cv_df$RMSE[irow] <- sqrt(mean((obs_cpue - pred_cpue)^2)) 
    
    #Update progress
    print(paste0("Done with: ", cv_df$species[irow], ", ", 
                 ifelse(cv_df$depth[irow], "Depth, ", "No Depth, "),
                 "Fold Number ", cv_df$fold[irow]))
  }
}

##################################################
####   Calculate how many of the folds converged with < 1E-5
##################################################
tidyr::spread(data = aggregate(max_grad ~ species + depth,
                               data = cv_df,
                               FUN = function(x) sum(x < 1E-4)),
              key = "depth",
              value = "max_grad")

##################################################
####  Calculate mean out-of-bag predicted NLL across folds
##################################################
prednll <- tidyr::spread(data = aggregate(pred_nll ~ species + depth,
                                          data = cv_df,
                                          FUN = function(x) round(mean(x)),
                                          subset = max_grad < 1E-4),
                         key = "depth",
                         value = "pred_nll")
prednll$depth_in_model <- apply(X = prednll[, c("FALSE", "TRUE")],
                                MARGIN = 1,
                                FUN = function(x) c("FALSE", "TRUE")[which.min(x)] )

RMSE <- tidyr::spread(data = aggregate(RMSE ~ species + depth,
                                       data = cv_df,
                                       FUN = function(x) round(sqrt(mean(x^2))),
                                       subset = max_grad < 1E-4),
                      key = "depth",
                      value = "RMSE")

##################################################
####   Create the result object that would go into the optimizations
##################################################
N <- nrow(Extrapolation_depths)
D_gct = Index <- array(dim = c(N, ns, n_years), 
                       dimnames = list(NULL, which_spp, NULL))

for(ispp in 1:ns){
  depth_in_model <- prednll$depth_in_model[ispp]
  result_dir = paste0(VAST_dir, prednll$species[ispp],
                      ifelse(depth_in_model, "_depth", ""),
                      "/")
  
  load(paste0(result_dir, "/fit.RData"))
  
  D_gct[,ispp,] = fit$Report$D_gct[, 1 , which_years]
  Index[,ispp,] = fit$Report$Index_gctl[, 1, which_years, 1]
  
  print(result_dir)
}

##################################################
####   Save
##################################################
save("RMSE", file = paste0(github_dir, "/data/rmse_VAST_models.RData"))
save("prednll", file = paste0(github_dir, "/data/prednll_VAST_models.RData"))
save("D_gct", file = paste0(github_dir, "/data/fit_density.RData") )
save("Index", file = paste0(github_dir, "/data/fit_index.RData") )

