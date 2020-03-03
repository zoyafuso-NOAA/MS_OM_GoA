###############################
## Model Settings for VAST
###############################
library(VAST)

# Version of VAST
rm(list = ls())

setwd( 'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/')

modelno = '5c'
if(!dir.exists(paste0(getwd(), '/VAST_output', modelno, '/'))) {
  dir.create(paste0(getwd(), '/VAST_output', modelno, '/'))
}

Version = get_latest_version( package="VAST" )

## Model settings
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) 
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,0)   


##Derived Products: We also decide on which post-hoc calculations to include in the output
Options =  c("SD_site_density"=1, 
             "SD_site_logdensity"=1, 
             "Project_factors" = 1)

save.image(paste0(getwd(), '/VAST_output', modelno, '/Model_Settings.RData') )
