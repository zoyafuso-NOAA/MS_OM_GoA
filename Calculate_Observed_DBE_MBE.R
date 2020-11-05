###############################################################################
## Project:     Calculate Observed CVs
## Author:      Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description: Calculate the observed CV in the surveys by boat effort
###############################################################################
rm(list = ls())

##################################################
####  Import libraries
##################################################
library(tidyr)

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/MS_OM_GoA/data/")

VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/Single_Species/"

##################################################
####  Constants
##################################################
sci_names <- c("Atheresthes stomias", "Gadus chalcogrammus",
               "Gadus macrocephalus", "Glyptocephalus zachirus" , 
               "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
               "Lepidopsetta bilineata", "Lepidopsetta polyxystra", 
               "Microstomus pacificus", "Sebastes alutus", 
               "Sebastes brevispinis", "Sebastes polyspinis", 
               "Sebastes variabilis", "Sebastolobus alascanus")

which_years <- c(1996, 1999, seq(from = 2003, to = 2019, by = 2))

# haul_n <- unique(DBE_CV[, c("YEAR", "HAUL_COUNT")])
# haul_n <- haul_n[order(haul_n$YEAR), ]
# three_boat_yrs <- haul_n$YEAR[haul_n$HAUL_COUNT > 700]
# two_boat_yrs <- haul_n$YEAR[haul_n$HAUL_COUNT < 700]

# Year_Set <- 1996:2019
# Years2Include <- c(1,  4,  8, 10, 12, 14, 16, 18, 20, 22, 24)

##################################################
####  Import data
##################################################
DBE_CV <- readRDS(paste0(github_dir, 'GOA_biomass_indices_wnames.rds'))

##################################################
####  Subset DBE data to species and years of interest
##################################################
DBE_CV <- subset(x = DBE_CV, 
                 subset = SPECIES_NAME %in% sci_names & YEAR %in% which_years,
                 select = c("YEAR", "SPECIES_NAME", 
                            "MEAN_WGT_CPUE", "VAR_WGT_CPUE"))

DBE_CV$CV_WGT_CPUE <- sqrt(DBE_CV$VAR_WGT_CPUE) / DBE_CV$MEAN_WGT_CPUE
DBE_CV <- DBE_CV[, c(1,2,5)]

DBE_CV <- tidyr::spread(data = DBE_CV, 
                        value = CV_WGT_CPUE, 
                        key = SPECIES_NAME)


VAST_CV <- VAST_CV_depth <- data.frame()

for (spp in sci_names) {
  for (depth_in_model in c(F, T)) 
  {
    filename <- paste0(VAST_dir, 
                       spp,
                       ifelse(depth_in_model, "_depth", ""),
                       "/diagnostics/Index/Table_for_SS3.csv")
    
    if (file.exists(filename)) {
      Index <- read.csv(filename)
      Index <- subset(Index, 
                      subset = Year %in% which_years,
                      select = c("Year", "SD_log"))
      Index$spp <- spp
      
      var_name <- paste0("VAST_CV", ifelse(depth_in_model, "_depth", "") )
      
      assign(x = var_name, value = rbind(get(var_name), Index))
      
    }
  }
}

VAST_CV <- tidyr::spread(VAST_CV, value = "SD_log", key = spp)
VAST_CV_depth <- tidyr::spread(VAST_CV_depth, value = "SD_log", key = spp)
names(DBE_CV)[1] <- "Year"

master_CV <- rbind(cbind(Method = "DBE", DBE_CV), 
                   cbind(Method = "MBE", VAST_CV),
                   cbind(Method = "MBE_depth", VAST_CV_depth))

names(master_CV)[-(1:2)] <- gsub(x = sci_names,
                                 pattern = " ",
                                 replacement = "_") 

par(mar = c(0,3,0,1),
    oma = c(2,0, 2, 0),
    mfrow = c(5,3))

for (spp in sci_names) {
  
  spp_name = gsub(x = spp,
                  pattern = " ",
                  replacement = "_")
  
  ylim_ <- max(master_CV[, spp_name]) * 1.1
  
  plot(1, 
       type = "n",
       xlim = c(0.5, 3.5),
       ylim = c(0, ylim_),
       axes = F,
       ann = F)
  axis(side = 2, las = 1)
  box()
  mtext(side = 3, 
        line = -2,
        text = spp,
        font = 3)

  
  boxplot(formula(paste(spp_name, " ~ Method")),
          data = master_CV[, c("Method", "Year", spp_name)],
          main = spp,
          add = T,
          axes = F,
          pch = 16,
          cex = 2)
  
  if(spp %in% sci_names[length(sci_names):(length(sci_names)-2)]) 
    axis(side = 1, at = 1:3, labels = c("DBE", "VAST", "VAST_depth") )
}

