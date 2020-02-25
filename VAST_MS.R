########################
## Build and Run Multispecies VAST
########################

rm(list = ls())

# Getting started
library(TMB)               # Can instead load library(TMBdebug)
library(VAST)

modelno = "4b"

setwd(paste0('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/VAST_output', modelno))

## Import Settings
load('Model_Settings.RData')
load('Spatial_Settings.RData')

# Build and run model

## Build model: To estimate parameters, we first build a list of data-inputs used for parameter estimation.  `make_data` has some simple checks for buggy inputs, but also please read the help file `?make_data`.  
#We then build the TMB object.

TmbData = make_data("Version"=Version, 
                    "FieldConfig"=FieldConfig, 
                    "OverdispersionConfig"=OverdispersionConfig, 
                    "RhoConfig"=RhoConfig, 
                    "ObsModel_ez" = c(PosDist = 2, Link = 0), 
                    "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                    "b_i"=Data_Geostat[,'Catch_KG'], 
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
                    "s_i"=Data_Geostat[,'knot_i']-1, 
                    "t_i"=Data_Geostat[,'Year'], 
                    "spatial_list"=Spatial_List, 
                    "Options"=Options ,
                    formula = "Catch_KG ~ DEPTH + DEPTH2",
                    covariate_data = cbind(Data_Geostat[,c('Lat', 'Lon', 
                                                           'DEPTH', 'DEPTH2',
                                                           'Catch_KG')], 
                                           Year = NA)
)


TmbList = make_model("TmbData"=TmbData, 
                     "RunDir"= getwd(), 
                     "Version"=Version, 
                     "RhoConfig"=RhoConfig, 
                     "loc_x"=Spatial_List$loc_x, 
                     "Method"=Spatial_List$Method)
Obj = TmbList[["Obj"]]

## Estimate fixed effects and predict random effects: Next, we use a gradient-based nonlinear minimizer to identify maximum likelihood estimates for fixed-effects

Opt = TMBhelper::fit_tmb( obj=Obj, 
                          lower=TmbList[["Lower"]], 
                          upper=TmbList[["Upper"]], 
                          getsd=TRUE,
                          getJointPrecision = TRUE,
                          savedir=getwd(), 
                          bias.correct=F, 
                          quiet = T,
                          bias.correct.control=list(
                            sd=F, split=NULL, 
                            nsplit=1, vars_to_correct="Index_cyl"), 
                          newtonsteps=1 )


#Finally, we bundle and save output
Report = Obj$report()

Save = list('Obj' = Obj,"Opt"=Opt, "Report"=Report, "TmbData"=TmbData, 
            'Spp' = unique(Data_Geostat$spp), "ParHat"=Obj$env$parList(Opt$par),
            'TmbList' = TmbList)

save(Save, file="VAST_MS_GoA_Run.RData")
