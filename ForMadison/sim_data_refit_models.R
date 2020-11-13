###############################################################################
## Project:      
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: 
## Description:  
###############################################################################
rm(list = ls())

##################################################
####   Load packages 
##################################################
library(VAST)

##################################################
####   Set up model settings
####   Create new directories for each model
####   Set up directories
##################################################
github_dir <- "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/"
github_dir2 <- "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/"
VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/"

##################################################
####   Load data
##################################################
master_data <- read.csv(file = paste0(github_dir, "data/GOA_multspp.csv") )
load(paste0(github_dir2, "data/RMSE_VAST_models.RData"))
load(paste0(github_dir2, "data/optimization_data.RData"))
load(paste0(VAST_dir, "sim_data.RData"))
source(paste0(github_dir, "fit_model_X_GTP.R"))

##################################################
####   Loop over species
##################################################

for (ispp in 1:14) {
  
  depth_in_model <- RMSE[ispp, "depth_in_model"]
  
  #Temporary directory
  spp_dir <- paste0(VAST_dir, "Single_Species/", sci_names[ispp],
                    ifelse(depth_in_model, "_depth", ""), "/" )
  
  ##################################################
  ####   Subset species
  ##################################################
  data <- subset(master_data, 
                 SPECIES_NAME == sci_names[ispp])
  
  ##################################################
  ####   Spatial settings: The following settings define the spatial resolution 
  ####   for the model, and whether to use a grid or mesh approximation
  ####   Stratification for results
  ##################################################
  Method <- "Mesh"
  strata.limits <- data.frame("STRATA" = c("All_areas"),
                              "west_border" = -Inf,
                              "east_border" = Inf)
  
  settings <- FishStatsUtils::make_settings( 
    n_x = 500,   # Number of knots
    Region = "User",
    purpose = "index2",
    strata.limits = strata.limits, 
    bias.correct = FALSE,
    FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, 
                    "Omega2" = 1, "Epsilon2" = 1),
    RhoConfig = c("Beta1" = 0, "Beta2" = 0, 
                  "Epsilon1" = 0, "Epsilon2" = 0),
    OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0), 
    "Options" = c("Calculate_Range" = F, 
                  "Calculate_effective_area" = F),
    ObsModel = c(2, 0),
    use_anisotropy = T)
  
  ##################################################
  ####   Import "true" and not interpolated covariate 
  ####   data if using depth covariates
  ##################################################
  
  load( paste0(github_dir, "data/Extrapolation_depths.RData"))
  
  n_g <- nrow(Extrapolation_depths) #number of grid cells
  n_t <- diff(range(Data_Geostat$Year)) + 1 #Number of total years
  n_p <- 2 #two density covariates
  
  X_gtp <- array(dim = c(n_g, 1, n_t, n_p) )
  for (i in 1:n_t) {
    X_gtp[,,i,] <- 
      as.matrix(Extrapolation_depths[,c("LOG_DEPTH_EFH_CEN", 
                                        "LOG_DEPTH_EFH_CEN_SQ")])
  }
  
  for (iter in 1:10) {
    
    ##################################################
    ####   Create result directory
    ####   copy vast version cpp, dll, and o files
    ##################################################
    output_dir <- paste0(spp_dir, "sim_runs/iter_", iter)
    if (!dir.exists(output_dir)) dir.create(path = output_dir, recursive = T)
    
    file.copy(from = dir(spp_dir, pattern = "VAST", full.names = T),
              to = output_dir)
    
    ##################################################
    ####   Prepare the dataframe for catch-rate data in the VAST format
    ####   Replace "Catch_KG" with the simulated data
    ##################################################
    Data_Geostat <- data.frame( "spp" = data$SPECIES_NAME,
                                "Year" = data$YEAR,
                                "Catch_KG" = sim_data[iter, , ispp],
                                "AreaSwept_km2" = data$EFFORT,
                                "Vessel" = 0,
                                "Lat" = data$LATITUDE,
                                "Lon" = data$LONGITUDE, 
                                stringsAsFactors = T)
    
    if (depth_in_model){
      Data_Geostat[, c("LOG_DEPTH", 
                       "LOG_DEPTH2") ] = data[, c("LOG_DEPTH_EFH_CEN", 
                                                  "LOG_DEPTH_EFH_CEN_SQ")]
    }
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    if (!depth_in_model) {
      fit = FishStatsUtils::fit_model( 
        "settings" = settings,
        "working_dir" = output_dir,
        "Lat_i" = Data_Geostat[, "Lat"],
        "Lon_i" = Data_Geostat[, "Lon"],
        "t_i" = Data_Geostat[, "Year"],
        "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
        "b_i" = Data_Geostat[, "Catch_KG"],
        "a_i" = Data_Geostat[, "AreaSwept_km2"],
        "v_i" = Data_Geostat[, "Vessel"],
        "max_cells" = Inf,
        "getJointPrecision" = TRUE,
        "newtonsteps" = 1,
        "test_fit" = F,
        "Options" = c("Calculate_Range" = F, 
                      "Calculate_effective_area" = F),
        "input_grid" = Extrapolation_depths)
    }
    
    if (depth_in_model) {
      fit = fit_model( 
        "settings" = settings,
        "working_dir" = output_dir,
        "Lat_i" = Data_Geostat[, "Lat"],
        "Lon_i" = Data_Geostat[, "Lon"],
        "t_i" = Data_Geostat[, "Year"],
        "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
        "b_i" = Data_Geostat[, "Catch_KG"],
        "a_i" = Data_Geostat[, "AreaSwept_km2"],
        "v_i" = Data_Geostat[, "Vessel"],
        "max_cells" = Inf,
        "getJointPrecision" = TRUE,
        "newtonsteps" = 1,
        "test_fit" = F,
        ##Additional arguments for covariates
        "X1_formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
        "X2_formula" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
        "covariate_data" = cbind(Data_Geostat[,c("Lat", 
                                                 "Lon",
                                                 "LOG_DEPTH",
                                                 "LOG_DEPTH2",
                                                 "Catch_KG")],
                                 Year = NA),
        "X_gtp" = X_gtp,
        "input_grid" = Extrapolation_depths
      )
    }
    
    ##################################################
    ####   Save model
    ##################################################
    save(list = c("fit"), 
         file = paste0(output_dir, "/fit.RData"))
  }
}
