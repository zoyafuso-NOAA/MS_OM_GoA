########################
## Build and Run Multispecies VAST
########################

rm(list = ls())

# Getting started
#utils::install.packages( "https://inla.r-inla-download.org/R/stable/bin/windows/contrib/3.5/INLA_18.07.12.zip" )
#library(INLA)
library(TMBdebug)
#library(devtools)
#devtools::install_local("C:/Users/Zack Oyafuso/Downloads/FishStatsUtils-2.5.0")
library(VAST)

modelno = "6c"

setwd(paste0('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/VAST_output', modelno))
# setwd(paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/VAST_output', modelno))

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

#Add "true" and not interpolated covariate data
load('../Extrapolation_depths.RData')
X_gtp = array(dim = c(TmbData$n_g,TmbData$n_t, TmbData$n_p) )
for(i in 1:TmbData$n_t) {
  X_gtp[,i,] = as.matrix(Extrapolation_depths[,c('DEPTH', 'DEPTH2')])
}

TmbData$X_gtp = X_gtp

TmbList = make_model("TmbData"=TmbData, 
                     "RunDir"= getwd(), 
                     "Version"=Version, 
                     "RhoConfig"=RhoConfig, 
                     "loc_x"=Spatial_List$loc_x, 
                     "Method"=Spatial_List$Method)
Obj = TmbList[["Obj"]]

## Estimate fixed effects and predict random effects: Next, we use a gradient-based nonlinear minimizer to identify maximum likelihood estimates for fixed-effects
bias_correct = F
Opt = TMBhelper::fit_tmb( obj=Obj, 
                          lower=TmbList[["Lower"]], 
                          upper=TmbList[["Upper"]], 
                          getsd=TRUE,
                          getJointPrecision = FALSE,
                          savedir=getwd(), 
                          bias.correct=bias_correct, 
                          quiet = T, 
                          newtonsteps=1 )


#Finally, we bundle and save output
Report = Obj$report()

Save = list('Obj' = Obj,"Opt"=Opt, "Report"=Report, "TmbData"=TmbData, 
            'Spp' = levels(Data_Geostat$spp), "ParHat"=Obj$env$parList(Opt$par),
            'TmbList' = TmbList)

save(Save, file="VAST_MS_GoA_Run.RData")
