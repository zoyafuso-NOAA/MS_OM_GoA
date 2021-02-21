###############################################################################
## Project:       VAST diagnostics plots
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate various diagnostic plots
##                Data and knot locations
##                Encounter probability
##                QQ plots
##                Pearson Residuals 
##                Mean density andyear-specific density, 
##                Spatial (omega) and spatiotemporal (epsilon) effects, 
##                Index of abundance
##                Covariance
##                Factors loadings
##                Anisotropy
###############################################################################

##################################################
#### Set up directories     
##################################################
rm(list = ls())
which_machine <- c("Zack_PC" = 1, "Zack_GI_PC" = 2)[2]

github_dir <- c("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/")[which_machine]
VAST_dir <- c("C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Single_Species/",
              "G:/Oyafuso/VAST_Runs_EFH/Single_Species/")[which_machine]

##################################################
####  Import Libraries
##################################################
library(VAST)
library(raster)
library(sp)
library(rgdal)
library(RColorBrewer)
library(plotrix)
library(rnaturalearth)

##################################################
####  Constants
##################################################
sci_names <- c("Sebastes polyspinis", "Sebastes variabilis", 
               "Sebastes brevispinis", "Microstomus pacificus",
               "Lepidopsetta polyxystra", "Lepidopsetta bilineata",
               "Hippoglossus stenolepis", "Hippoglossoides elassodon",
               "Glyptocephalus zachirus", "Gadus macrocephalus",
               "Gadus chalcogrammus", "Sebastes B_R", "Sebastes alutus",
               "Atheresthes stomias", "Sebastolobus alascanus",
               "Anoplopoma fimbria", "Beringraja spp.", "Octopus spp.",
               "Pleurogrammus monopterygius", "Sebastes borealis",
               # "Sebastes ruberrimus", 
               "Sebastes variegatus", "Squalus suckleyi")

sci_names_filename <- gsub(x = sci_names, 
                           pattern = "\\.",
                           replacement = "")

for (depth_in_model in c(T, F)) {
  for (ispp in sci_names_filename[1]) {
    
    result_dir <- paste0(VAST_dir, ispp, 
                         ifelse(depth_in_model,  "_depth", ""), "/")
    
    fun_dir <- paste0(github_dir, "diagnostics/")
    
    if(!dir.exists(paste0(result_dir, "diagnostics/")) ) 
      dir.create(paste0(result_dir, "diagnostics/"))
    
    ##################################################
    ####  Load VAST fit and import customized plot functions 
    ##################################################
    load(paste0(result_dir, "/fit.RData"))
    
    ##################################################
    #### Extract objects from fitted object   
    #### Set up constants
    ##################################################
    data_geostat <- fit$data_frame
    names(data_geostat)[c(1:2, 5:7)] <- c("Lat", "Lon", 
                                          "Catch_KG","Year", "spp")
    
    opt <- fit$parameter_estimates
    report <- fit$Report
    tmbdata <- fit$data_list
    
    year_set <- seq(min(data_geostat[, "Year"]),
                    max(data_geostat[, "Year"]))
    years_included <- which( year_set %in% sort(unique(data_geostat[, "Year"])))
    
    index_ests <- report$Index_cyl[, years_included, 1]
    which_sds <- attributes(opt$SD$value)$names == "Index_ctl"
    index_sds <- matrix(data = opt$SD$sd[which_sds], 
                        nrow = 1 )[, years_included]
    
    ############################################
    ## Index of abundance
    ## Compare VAST estimates with Design-Based (Stratified RS) estimates
    ############################################
    Index = plot_biomass_index( DirName = paste0(result_dir, "diagnostics/"), 
                                TmbData = tmbdata, 
                                Sdreport = opt[["SD"]], 
                                Year_Set = year_set, 
                                Years2Include = years_included, 
                                strata_names = strata.limits <- data.frame(
                                  "STRATA" = c("All_areas"),#, "west_of_140W"),
                                  "west_border" = c(-Inf),#, -Inf),
                                  "east_border" = c(Inf)#, -140)
                                )[,1], 
                                use_biascorr = TRUE, 
                                category_names = levels(data_geostat[, "spp"]))
    
    ############################################
    ## Direction of "geometric anisotropy"
    ## We can visualize which direction has faster or slower decorrelation 
    ## (termed "geometric anisotropy")
    ############################################
    aniso <- plot_anisotropy( FileName = paste0(result_dir, 
                                                "diagnostics/Aniso.png"), 
                              Report = report, 
                              TmbData = tmbdata )
  }
  
  
}
