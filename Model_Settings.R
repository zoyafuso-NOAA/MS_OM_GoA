###############################
## Model Settings for VAST
###############################

# Version of VAST
rm(list = ls())

modelno = '3a'
Version = get_latest_version( package="VAST" )

## Model settings
FieldConfig = c("Omega1"=3, "Epsilon1"=3, "Omega2"=3, "Epsilon2"=3) 
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,0)   

##Derived Products: We also decide on which post-hoc calculations to include in the output
Options =  c("SD_site_density"=1, 
             "SD_site_logdensity"=0, 
             "Calculate_Range"=0, 
             "Calculate_evenness"=0, 
             "Calculate_effective_area"=0, 
             "Calculate_Cov_SE"=0, 
             'Calculate_Synchrony'=0, 
             'Calculate_Coherence'=0)

save.image('Model_Settings.RData')
