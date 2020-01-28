
rm(list=ls())  #clears history
gc()


####### Setting lists of required packages & installing it
list_rpackage <- c("raster", "rgdal", "maptools", "gstat", "rgeos", "proj4", 
                   "PBSmapping", "sp", "ncdf4", "INLA", "gridExtra",
                   "splancs", "automap", "gtools", "reshape", "fields", 
                   "dplyr", "ggplot2", "reshape2", "tidyverse", 'Sumfish')
which_not_installed <- which(list_rpackage %in% rownames(installed.packages()) == FALSE)
if(length(which_not_installed)>1)
{
  install.packages(list_rpackage[which_not_installed],dep=TRUE)
  if ("INLA" %in% list_rpackage[which_not_installed]) install.packages("INLA", repos=c(getOption("repos"), 
                                                                                       INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

# loading the packages
lapply(list_rpackage, require, character.only = TRUE)  

#run this once
 original_dataset <- sumHaul(getRacebase(year=c(1982, 2019), 
                                         survey = 'EBS_SHELF'))
 write.csv(original_dataset, 'original_dataset.csv')
 All_data <- original_dataset  #making safe copy of data


df = data.frame()

for(spp in c(10110, 10130, 10210)){
  temp_df = read.csv(paste0(spp, 'final.csv'))
  df = rbind(df,  data.frame(cbind(spp = spp, temp_df)))
}
rm(temp_df, spp)

for(id in unique(df$STATIONID) ){
  temp_df = subset(df, STATIONID == id)
}
