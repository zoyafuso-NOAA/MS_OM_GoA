###############################################################################
## Project:       Knit together VAST results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate CVs for VAST Runs, save to file 
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
github_dir <- "C:/Users/Zack Oyafuso/Documents/GitHub/Optimal_Allocation_GoA/"
VAST_dir <- "C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Single_Species/"

##################################################
####    Import required packages
##################################################
library(FishStatsUtils)

##################################################
####    Load Data
##################################################
load(paste0(github_dir, "data/RMSE_VAST_models.RData"))
load(paste0(github_dir, "data/optimization_data.RData"))

##################################################
####    Result Object
##################################################
vast_index <- data.frame()

for (ispp in 1:ns_all) {
  #Load vast output
  temp_index <- read.csv(
    paste0(VAST_dir, 
           common_names_all[ispp], 
           ifelse(RMSE$depth_in_model[ispp], "_depth", ""),
           "/diagnostics/Table_for_SS3.csv"))[years_included,]
  
  #Record CV 
  vast_index = rbind(vast_index,
                     with( temp_index, 
                           data.frame(
                             spp = common_names_all[ispp],
                             year = Year,
                             est = Estimate_metric_tons, 
                             cv = SD_mt / Estimate_metric_tons)
                     )
  )
}

##################################################
####    Plot for show
##################################################
par(mar = c(3,12,1,1))
boxplot(cv ~ spp, 
        data = vast_index,
        horizontal = T,
        las = 1,
        ann = F,
        ylim = c(0, 0.5))

##################################################
####    Save
##################################################
save(list = "vast_index", 
     file = paste0(github_dir, "data/vast_index.RData"))
