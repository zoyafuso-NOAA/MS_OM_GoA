###############################################################################
## Project:      Univariate VAST model runs
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
##               Jim Thorson"s VAST wiki example
##           (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:  Run single-species VAST models wit and without depth
##               as a covariate. Run 10-fold Cross Validation for each Model
##
##               Software versions:
##               R version 4.0.2 (2020-06-22)
##               VAST version 3.6.1
##               FishStatsUtils 2.8.0
##               VAST_v12_0_0.cpp
###############################################################################
rm(list = ls())

##################################################
####   Since the analysis is not contained within a directory, path names
####   are set up here to be used throughout. The script is housed in a a 
####   github directory whereas the model output is sent to Z. Oyafuso's G
####   Drive due to storage limits on github.
##################################################
which_machine <- c("Zack_PC" = 1, "Zack_GI_PC" = 2)[1]

github_dir <- c("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/")[which_machine]
VAST_dir <- c("C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Single_Species/",
              "G:/Oyafuso/VAST_Runs_EFH/Single_Species/")[which_machine]

if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)

##################################################
####   Load packages, make sure versions are consistent
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


##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
master_data <- read.csv(file = paste0(github_dir, "data/GOA_multspp.csv") )

#################################################
## Loop over species to fit models with and without depth covariates
#################################################
spp_names <- sort(unique(master_data$COMMON_NAME))

for (ispp in spp_names[1]) {
  for (depth_in_model in c(F, T)[2]) {
    
    ##################################################
    ## Create directory to store model results
    ##################################################
    result_dir <- paste0(VAST_dir, 
                         ispp, 
                         ifelse(test = depth_in_model,  
                                yes = "_depth", 
                                no = ""), 
                         "/")
    
    if (!dir.exists(result_dir)) dir.create(result_dir)  
    
    ##################################################
    ## Create diagnostics objects
    ##################################################
    if(!dir.exists(paste0(result_dir, "/diagnostics"))) {
      dir.create(paste0(result_dir, "/diagnostics"))
    }
    
    ##################################################
    ####   Subset species
    ##################################################
    data <- subset(master_data, 
                   COMMON_NAME == ispp)
    
    ##################################################
    ####   Prepare the dataframe for catch-rate data in the VAST format
    ##################################################
    Data_Geostat <- data.frame( spp = data$SPECIES_NAME,
                                Year = data$YEAR,
                                Catch_KG = data$WEIGHT,
                                AreaSwept_km2 = data$EFFORT,
                                Lat = data$LATITUDE,
                                Lon = data$LONGITUDE, 
                                stringsAsFactors = T)
    
    Data_Geostat[, c("LOG_DEPTH", "LOG_DEPTH2") ] <-
      data[, c("LOG_DEPTH_EFH_CEN", "LOG_DEPTH_EFH_CEN_SQ")]
    
    ##################################################
    ####   Assign 10 fold partitions of the data
    ##################################################
    n_fold <- 10
    years <- paste0(unique(Data_Geostat$Year))
    NTime <- length(unique(Data_Geostat$Year))
    
    #Create unique stationID from the latlon. To make sure the ids are unique,
    #we use the table function to make sure there are 7900 records (as of 2019)=
    Data_Geostat$latlon <- paste0(Data_Geostat$Lat, Data_Geostat$Lon)
    table(table(Data_Geostat$latlon))
    
    #split Data_Geostat by year, then on each year-split, randomly assign 
    #fold numbers to the each unique station
    set.seed(2342)
    foldno <- lapply(
      #Split Data_Geostat by Year
      X = split.data.frame(Data_Geostat, 
                           f = Data_Geostat$Year),
      
      #For each year split, randomly assign fold numbers so that each year is 
      #equally split into n_folds folds
      FUN = function(x) {
        unique_loc <- unique(x$latlon)
        fold_no <- sample(x = 1:n_fold, 
                          size = length(unique_loc), 
                          replace = T)
        return(split(unique_loc, fold_no))
      })
    
    #Attach fold number to the Data_Geostat
    for (iyear in years) {
      for (ifold in paste(1:n_fold)) {
        Data_Geostat[Data_Geostat$latlon %in% foldno[[iyear]][[ifold]] , 
                     "fold"] = as.integer(ifold) 
      }
    }
    
    #Columns should roughly have the same number of samples
    # table(Data_Geostat$fold, Data_Geostat$Year)
    
    ##################################################
    ####   Spatial settings: The following settings define the spatial resolution 
    ####   for the model, and whether to use a grid or mesh approximation
    ####   Stratification for results
    ##################################################
    settings <- FishStatsUtils::make_settings( 
      Version = "VAST_v12_0_0",
      n_x = 500,   # Number of knots
      Region = "User", #User inputted extrapolation grid
      purpose = "index2",
      fine_scale = TRUE,
      strata.limits =  data.frame("STRATA" = c("All_areas"),
                                  "west_border" = -Inf,
                                  "east_border" = Inf), 
      bias.correct = FALSE,
      FieldConfig = c(
        "Omega1" = 1,   #Spatial random effect on occurence 
        "Epsilon1" = 1, #Spatiotemporal random effect on occurence 
        "Omega2" = 1,   #Spatial random effect on positive response 
        "Epsilon2" = 1  #Spatiotemporal random effect on positive response
      ), 
      RhoConfig = c("Beta1" = 0, 
                    "Beta2" = 0, 
                    "Epsilon1" = 0, 
                    "Epsilon2" = 0), #Each year is a fixed effect
      OverdispersionConfig = c("Eta1" = 0, 
                               "Eta2" = 0), #Turn off overdispersion 
      "Options" = c("Calculate_Range" = F, 
                    "Calculate_effective_area" = F),
      
      ObsModel = c(2, 1),
      max_cells = Inf,
      use_anisotropy = T)
    
    ##################################################
    ####   Import "true" and not interpolated covariate 
    ####   data if using depth covariates
    ##################################################
    load( paste0(github_dir, "data/Extrapolation_depths.RData"))
    
    n_g <- nrow(Extrapolation_depths) #number of grid cells
    n_t <- diff(range(Data_Geostat$Year)) + 1 #Number of total years
    n_p <- 2 #two density covariates
    
    X_gtp <- array(dim = c(n_g, n_t, n_p) )
    for (i in 1:n_t) {
      X_gtp[, i, ] <- 
        as.matrix(Extrapolation_depths[,c("LOG_DEPTH_EFH_CEN", 
                                          "LOG_DEPTH_EFH_CEN_SQ")])
    }
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    fit = FishStatsUtils::fit_model( 
      "settings" = settings,
      "working_dir" = result_dir,
      "Lat_i" = Data_Geostat[, "Lat"],
      "Lon_i" = Data_Geostat[, "Lon"],
      "t_i" = Data_Geostat[, "Year"],
      "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
      "b_i" = Data_Geostat[, "Catch_KG"],
      "a_i" = Data_Geostat[, "AreaSwept_km2"],
      "getJointPrecision" = TRUE,
      "newtonsteps" = 1,
      "test_fit" = F,
      "input_grid" = Extrapolation_depths,
      
      ## Extra arguments if depth covariates are included
      "X1_formula" =  switch(depth_in_model,
                             "TRUE" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                             "FALSE" = "~0"),
      "X2_formula" =  switch(depth_in_model,
                             "TRUE" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                             "FALSE" = "~0"),
      "covariate_data" = switch(depth_in_model,
                                "TRUE" = cbind(Data_Geostat[,c("Lat",
                                                               "Lon",
                                                               "LOG_DEPTH",
                                                               "LOG_DEPTH2",
                                                               "Catch_KG")],
                                               Year = NA),
                                "FALSE" = NULL),
      "X_gtp" = switch(depth_in_model,
                       "TRUE" = X_gtp,
                       "FALSE" = NULL)
    )
    
    ## Diagnostics
    plot(x = fit,
         working_dir = paste0(result_dir, "diagnostics/"))
    
    # fit = switch(paste0(depth_in_model),
    #              "FALSE" = FishStatsUtils::fit_model( 
    #                "settings" = settings,
    #                "working_dir" = result_dir,
    #                "Lat_i" = Data_Geostat[, "Lat"],
    #                "Lon_i" = Data_Geostat[, "Lon"],
    #                "t_i" = Data_Geostat[, "Year"],
    #                "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
    #                "b_i" = Data_Geostat[, "Catch_KG"],
    #                "a_i" = Data_Geostat[, "AreaSwept_km2"],
    #                "getJointPrecision" = TRUE,
    #                "newtonsteps" = 1,
    #                "test_fit" = F,
    #                "input_grid" = Extrapolation_depths),
    #              
    #              "TRUE" = FishStatsUtils::fit_model( 
    #                "settings" = settings,
    #                "working_dir" = result_dir,
    #                "Lat_i" = Data_Geostat[, "Lat"],
    #                "Lon_i" = Data_Geostat[, "Lon"],
    #                "t_i" = Data_Geostat[, "Year"],
    #                "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
    #                "b_i" = Data_Geostat[, "Catch_KG"],
    #                "a_i" = Data_Geostat[, "AreaSwept_km2"],
    #                "getJointPrecision" = TRUE,
    #                "newtonsteps" = 1,
    #                "test_fit" = F,
    #                "input_grid" = Extrapolation_depths,
    #                
    #                ##Additional arguments for covariates
    #                "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
    #                "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
    #                "covariate_data" = cbind(Data_Geostat[,c("Lat",
    #                                                         "Lon",
    #                                                         "LOG_DEPTH",
    #                                                         "LOG_DEPTH2",
    #                                                         "Catch_KG")],
    #                                         Year = NA),
    #                "X_gtp" = X_gtp )
    # )
    
    ##################################################
    ####   Save
    ##################################################
    save(list = c("fit", "Data_Geostat"), 
         file = paste0(result_dir, "/fit.RData"))
    
    ##################################################
    ####   10-fold Cross Validation
    ##################################################
    n_fold <- 10
    for (fI in 1:n_fold) { 
      if (!dir.exists(paste0(result_dir, "CV_", fI))) {
        dir.create(paste0(result_dir, "CV_", fI))
        
        file.copy(from = paste0(result_dir, get_latest_version(), 
                                c(".cpp", ".dll", ".o")),
                  to = paste0(result_dir, "CV_", fI, "/", 
                              get_latest_version(), 
                              c(".cpp", ".dll", ".o")))
        
      }
    } 
    
    # Loop through partitions, refitting each time with a different PredTF_i
    for (fI in 1:n_fold ) {
      PredTF_i <- ifelse( test = Data_Geostat$fold == fI, 
                          yes = TRUE, 
                          no = FALSE )
      
      fit_new <- FishStatsUtils::fit_model( 
        "settings" = settings,
        "working_dir" = paste0(result_dir, "CV_", fI, "/"),
        "Lat_i" = Data_Geostat[, "Lat"],
        "Lon_i" = Data_Geostat[, "Lon"],
        "t_i" = Data_Geostat[, "Year"],
        "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
        "b_i" = Data_Geostat[, "Catch_KG"],
        "a_i" = Data_Geostat[, "AreaSwept_km2"],
        "getJointPrecision" = TRUE,
        "newtonsteps" = 1,
        "test_fit" = F,
        "input_grid" = Extrapolation_depths,
        "Parameters" = fit$ParHat,
        
        ## For cross validation runs, specify which folds are used in the 
        ## fitting process and which folds are withheld
        "PredTF_i" = PredTF_i, 
        
        ## Extra arguments if depth covariates are included
        "X1_formula" =  switch(depth_in_model,
                               "TRUE" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                               "FALSE" = "~0"),
        "X2_formula" =  switch(depth_in_model,
                               "TRUE" = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                               "FALSE" = "~0"),
        "covariate_data" = switch(depth_in_model,
                                  "TRUE" = cbind(Data_Geostat[,c("Lat",
                                                                 "Lon",
                                                                 "LOG_DEPTH",
                                                                 "LOG_DEPTH2",
                                                                 "Catch_KG")],
                                                 Year = NA),
                                  "FALSE" = NULL),
        "X_gtp" = switch(depth_in_model,
                         "TRUE" = X_gtp,
                         "FALSE" = NULL)
      )
      
      # fit_new = switch(paste0(depth_in_model),
      #                  "FALSE" = FishStatsUtils::fit_model( 
      #                    "settings" = settings,
      #                    "working_dir" = paste0(result_dir, "CV_", fI, "/"),
      #                    "Lat_i" = Data_Geostat[, "Lat"],
      #                    "Lon_i" = Data_Geostat[, "Lon"],
      #                    "t_i" = Data_Geostat[, "Year"],
      #                    "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
      #                    "b_i" = Data_Geostat[, "Catch_KG"],
      #                    "a_i" = Data_Geostat[, "AreaSwept_km2"],
      #                    "getJointPrecision" = TRUE,
      #                    "newtonsteps" = 1,
      #                    "test_fit" = F,
      #                    "input_grid" = Extrapolation_depths,
      #                    
      #                    "PredTF_i" = PredTF_i, 
      #                    "Parameters" = fit$ParHat,
      #                    "getsd" = T),
      #                  
      #                  "TRUE" = FishStatsUtils::fit_model( 
      #                    "settings" = settings,
      #                    "working_dir" = paste0(result_dir, "CV_", fI, "/"),
      #                    "Lat_i" = Data_Geostat[, "Lat"],
      #                    "Lon_i" = Data_Geostat[, "Lon"],
      #                    "t_i" = Data_Geostat[, "Year"],
      #                    "c_i" = as.numeric(Data_Geostat[, "spp"]) - 1,
      #                    "b_i" = Data_Geostat[, "Catch_KG"],
      #                    "a_i" = Data_Geostat[, "AreaSwept_km2"],
      #                    "getJointPrecision" = TRUE,
      #                    "newtonsteps" = 1,
      #                    "test_fit" = F,
      #                    "input_grid" = Extrapolation_depths,
      #                    
      #                    "PredTF_i" = PredTF_i, 
      #                    "Parameters" = fit$ParHat,
      #                    "getsd" = T,
      #                    
      #                    ##Additional arguments for covariates
      #                    "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
      #                    "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
      #                    "covariate_data" = cbind(Data_Geostat[,c("Lat",
      #                                                             "Lon",
      #                                                             "LOG_DEPTH",
      #                                                             "LOG_DEPTH2",
      #                                                             "Catch_KG")],
      #                                             Year = NA),
      #                    "X_gtp" = X_gtp )
      # )
      
      # Save fit 
      save(list = "fit_new",  
           file = paste0(result_dir, "CV_", fI, "/fit.RData"))
    }
  }
}
