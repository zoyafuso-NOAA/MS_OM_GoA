###############################################################################
## Project:      Simulate data from VAST
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
## Description:  Simulate densities according to the best-fitting VAST model
###############################################################################
rm(list = ls())

##################################################
####   Since the analysis is not contained within a directory, path names
####   are set up here to be used throughout. The script is housed in a a 
####   github directory whereas the model output is sent to Z. Oyafuso's G
####   Drive due to storage limits on github.
##################################################
which_machine <- c("Zack_PC" = 1, "Zack_GI_PC" = 2)[1]

github_dir <- paste0(c("C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/MS_OM_GoA/")
github_dir2 <- paste0(c("C:/Users/Zack Oyafuso/Documents/",
                        "C:/Users/zack.oyafuso/Work/")[which_machine], 
                      "GitHub/Optimal_Allocation_GoA/")
VAST_dir <- paste0(c("C:/Users/Zack Oyafuso/",
                     "C:/Users/zack.oyafuso/")[which_machine],
                   "Desktop/VAST_Runs/Simulate_Data/")

if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)

##################################################
####   Load packages, make sure versions are consistent
####   assign the version of the cpp file used
##################################################
library(VAST)

{
  switch(EXPR = R.version$version.string == "R version 4.0.2 (2020-06-22)",
         "TRUE" = print("R version is consistent with R version 4.0.2"),
         "FALSE" = print(
           paste0("R version ", 
                  R.version$version.string, 
                  " is not consistent with R version 4.0.2 (2020-06-22).",
                  "Update R version to be consistent")))
  
  switch(EXPR = packageVersion("VAST") == "3.6.1",
         "TRUE" = print("R package VAST is consistent with version 3.6.1"),
         "FALSE" = print(
           paste0("R version ", packageVersion("VAST"), 
                  " is not consistent with version 3.6.1.",
                  "Update R package VAST to be consistent")))
  
  switch(EXPR = packageVersion( "FishStatsUtils") == "2.8.0",
         "TRUE" = "R package FishStatsUtils is consistent with version 2.8.0",
         "FALSE" = print(
           paste0("R version ", packageVersion( "FishStatsUtils"), 
                  " is not consistent with version 2.8.0.",
                  "Update R package FishStatsUtils to be consistent")))
}

vast_cpp_version <- "VAST_v12_0_0"

##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
master_data <- read.csv(file = paste0(github_dir, "data/GOA_multspp.csv") )
load(paste0(github_dir2, "data/prednll_VAST_models.RData"))

##################################################
####   Import extrapolation grid with depth values
##################################################
load(paste0(github_dir, "data/Extrapolation_depths.RData"))

n_g <- nrow(Extrapolation_depths) #number of grid cells
n_t <- diff(range(master_data$YEAR)) + 1 #Number of total years
n_p <- 2 #two density covariates

X_gtp <- array(dim = c(n_g, n_t, n_p) )
for (i in 1:n_t) {
  X_gtp[, i, ] <- 
    as.matrix(Extrapolation_depths[,c("LOG_DEPTH_EFH_CEN", 
                                      "LOG_DEPTH_EFH_CEN_SQ")])
}

#################################################
## Loop over species to fit models with and without depth covariates
#################################################
for (irow in 1:nrow(prednll) ) {
  
  ##################################################
  ## Create directory to store model results
  ##################################################
  result_dir <- paste0(VAST_dir, prednll$species[irow], "/")
  if (!dir.exists(result_dir)) dir.create(result_dir)  
  
  ##################################################
  ####   Subset species
  ##################################################
  data <- subset(master_data, 
                 COMMON_NAME == prednll$species[irow])
  
  ##################################################
  #### Prediction Grid: df of the extrapolation grid to simulate data onto
  ##################################################
  extrapolation_grid_df <- data.frame()
  for(itime in unique(data$YEAR)) {
    extrapolation_grid_df <- 
      rbind(extrapolation_grid_df,
            data.frame(spp = prednll$species[irow],
                       Year = rep(itime, nrow(Extrapolation_depths)),
                       Catch_KG = mean(data$WEIGHT),
                       AreaSwept_km2 = Extrapolation_depths$Area_km2,
                       Lat = Extrapolation_depths$Lat,
                       Lon = Extrapolation_depths$Lon,
                       LOG_DEPTH = Extrapolation_depths$LOG_DEPTH_EFH_CEN,
                       LOG_DEPTH2 = Extrapolation_depths$LOG_DEPTH_EFH_CEN_SQ,
                       stringsAsFactors = T)
      )
  }
  
  ##################################################
  ####   Extract whether depth is included in model and field configurations
  ##################################################
  depth_in_model <- as.logical(prednll$depth_in_model[irow]) 
  
  load(paste0(dirname(VAST_dir), "/Single_Species/", prednll$species[irow], 
              ifelse(test = depth_in_model,
                     yes = "_depth",
                     no = ""),
              "/fit.RData"))
  original_fit <- fit
  rm(fit)
  
  temp_FieldConfig  <- original_fit$settings$FieldConfig
  
  
  ##################################################
  ####   Prepare the dataframe for catch-rate data in the VAST format
  ##################################################
  data_geostat <- data.frame( spp = data$COMMON_NAME,
                              Year = data$YEAR,
                              Catch_KG = data$WEIGHT,
                              AreaSwept_km2 = data$EFFORT,
                              Lat = data$LATITUDE,
                              Lon = data$LONGITUDE, 
                              stringsAsFactors = T)
  
  data_geostat[, c("LOG_DEPTH", "LOG_DEPTH2") ] <-
    data[, c("LOG_DEPTH_EFH_CEN", "LOG_DEPTH_EFH_CEN_SQ")]
  
  ###################################################
  ## Add New Points: set catch to NAs?
  ###################################################
  data_geostat_with_grid <- rbind(data_geostat,
                                  extrapolation_grid_df)
  
  ##################################################
  ####   Spatial settings: The following settings define the spatial resolution 
  ####   for the model, and whether to use a grid or mesh approximation
  ####   Stratification for results
  ##################################################
  settings <- FishStatsUtils::make_settings( 
    Version = vast_cpp_version,
    n_x = 500,   # Number of knots
    Region = "User", #User inputted extrapolation grid
    purpose = "index2",
    fine_scale = TRUE,
    strata.limits =  data.frame("STRATA" = c("All_areas"),
                                "west_border" = -Inf,
                                "east_border" = Inf), 
    bias.correct = FALSE,
    FieldConfig = temp_FieldConfig,
    "Options" = c("Calculate_Range" = F, 
                  "Calculate_effective_area" = F),
    
    ObsModel = c(2, 1),
    max_cells = Inf,
    use_anisotropy = T)
  
  ##################################################
  ####   Fit the model and save output
  ##################################################
  pred_TF <- rep(1, nrow(data_geostat_with_grid)); pred_TF[1:nrow(data)] <- 0
  
  fit = tryCatch( {
    switch(paste0(depth_in_model),
           "FALSE" = FishStatsUtils::fit_model( 
             "settings" = settings,
             "working_dir" = result_dir,
             "Lat_i" = data_geostat[, "Lat"],
             "Lon_i" = data_geostat[, "Lon"],
             "t_i" = data_geostat[, "Year"],
             "c_i" = as.numeric(data_geostat[, "spp"]) - 1,
             "b_i" = data_geostat[, "Catch_KG"],
             "a_i" = data_geostat[, "AreaSwept_km2"],
             "getJointPrecision" = TRUE,
             "newtonsteps" = 1,
             "test_fit" = F,
             "input_grid" = Extrapolation_depths, 
             "PredTF_i" = pred_TF,
             "Parameters" = original_fit$ParHat),
           
           "TRUE" = FishStatsUtils::fit_model( 
             "settings" = settings,
             "working_dir" = result_dir,
             "Lat_i" = data_geostat_with_grid[, "Lat"],
             "Lon_i" = data_geostat_with_grid[, "Lon"],
             "t_i" = data_geostat_with_grid[, "Year"],
             "c_i" = as.numeric(data_geostat_with_grid[, "spp"]) - 1,
             "b_i" = data_geostat_with_grid[, "Catch_KG"],
             "a_i" = data_geostat_with_grid[, "AreaSwept_km2"],
             "getJointPrecision" = TRUE,
             "newtonsteps" = 1,
             "test_fit" = F,
             "input_grid" = Extrapolation_depths[, c("Area_km2", "Lon", "Lat")],
             "Parameters" = original_fit$ParHat,
             
             ##Additional arguments for covariates
             "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
             "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
             "covariate_data" = cbind(data_geostat[, c("Lat",
                                                       "Lon",
                                                       "LOG_DEPTH",
                                                       "LOG_DEPTH2",
                                                       "Catch_KG")],
                                      Year = NA),
             "X_gtp" = X_gtp, 
             "PredTF_i" = pred_TF))},
    error = function(cond) return(NULL)
  )
  
  if( !is.null(fit) ) {
    save(fit, file = paste0(VAST_dir, prednll$species[irow], "/fit.RData"))
    
    sim_data <- array(data = NA, dim = c(nrow(Extrapolation_depths),
                                         length(unique(data$YEAR)),
                                         1000))
    
    if(!(vast_cpp_version %in% names(getLoadedDLLs()))) {
      dyn.load(paste0(result_dir, "/", vast_cpp_version, ".dll"))  
    }       
    for (isim in 1:1000) {
      Sim1 <- FishStatsUtils::simulate_data(fit = fit, 
                                            type = 1, 
                                            random_seed = isim)
      sim_data[, , isim] <- matrix(data = Sim1$b_i[pred_TF == 1] * 0.001, 
                                   nrow = nrow(Extrapolation_depths), 
                                   ncol = length(unique(data$YEAR)))
      if(isim%%100 == 0) print(paste0("Done with Iteration ", isim, ", ",
                                      prednll$species[irow]))
    }
    
    save(sim_data, file = paste0(result_dir, "/simulated_data.RData"))
    dyn.unload(paste0(result_dir, "/", vast_cpp_version, ".dll"))
  }
  
}

