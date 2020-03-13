####################
## Diagnostic plots
####################

rm(list = ls())

library(VAST); library(mvtnorm)

# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/diagnostics/')
setwd('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/')
source("summarize_covariance.R")
source("plot_residuals.R")
source("plot_factors.R")
source("plot_maps_density.R")
source("plot_variable_density.R")

# setwd('C:/Users/zack.oyafuso/Desktop/VAST_Runs')
setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs')
VAST_model = "6d"
load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj

DateFile = paste0('diagnostics/VAST_model', VAST_model, '/')
if(!dir.exists(DateFile)) dir.create(DateFile)

## Plot data: It is always good practice to conduct exploratory analysis of 
## data. Here, I visualize the spatial distribution of data. Spatio-temporal
## models involve the assumption that the probability of sampling a given
## location is statistically independent of the probability distribution for 
## the response at that location. So if sampling "follows" changes in density,
## then the model is probably not appropriate!

plot_data(Extrapolation_List=Extrapolation_List,
          Spatial_List=Spatial_List,
          Data_Geostat=Data_Geostat,
          PlotDir= DateFile)

## Convergence: Here I print the diagnostics generated during parameter 
## estimation,and I confirm that (1) no parameter is hitting an upper or lower
## bound and (2) the final gradient for each fixed-effect is close to zero. 
## For explanation of parameters, please see `?make_data`.
pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE','Upper',
                                         'final_gradient')] ) 


## Diagnostics for encounter-probability component: Next, we check whether observed 
## encounter frequencies for either low or high probability samples are within the 
## 95% predictive interval for predicted encounter probability
Enc_prob = plot_encounter_diagnostic( Report=Report, 
                                      cutpoints_z = seq(0, 1, length = 50),
                                      Data_Geostat=Data_Geostat, 
                                      DirName= DateFile)

## Diagnostics for positive-catch-rate component: We can visualize fit to residuals of 
## catch-rates given encounters using a Q-Q plot. A good Q-Q plot will have residuals 
## along the one-to-one line.

Q = plot_quantile_diagnostic( TmbData=TmbData, 
                              Report=Report, 
                              FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", 
                              FileName_Qhist="Q-Q_hist", 
                              DateFile=DateFile) 

## Diagnostics for plotting residuals on a map: Finally, we visualize residuals 
## on a map. To do so, we first define years to plot and generate plotting 
## inputs. useful plots by first determining which years to plot (`Years2Include`), and labels for each plotted 
## year (`Year_Set`)

# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"='Gulf_of_Alaska', 
                                 "spatial_list"=Spatial_List, 
                                 "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))


## We then plot Pearson residuals.  If there are visible patterns (areas 
## with consistently positive or negative residuals accross or within years) 
## then this is an indication of the model "overshrinking" results towards the 
## intercept, and model results should then be treated with caution.  
if(!dir.exists(paste0(DateFile, "Pearson_Residuals"))) 
  dir.create(paste0(DateFile, "Pearson_Residuals/"))

plot_residuals(Lat_i=Data_Geostat[,'Lat'], 
               Lon_i=Data_Geostat[,'Lon'], 
               TmbData=TmbData, 
               Report=Report, 
               Q=Q, 
               working_dir = paste0(DateFile, "Pearson_Residuals/"),
               xlim=MapDetails_List[["Xlim"]], 
               ylim=MapDetails_List[["Ylim"]], 
               Year_Set=Year_Set, Years2Include=Years2Include,                
		   zone=MapDetails_List[["Zone"]], 
               mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8,
               spatial_list = Spatial_List,
               extrapolation_list = Extrapolation_List)


## Plot spatial and spatio-temporal covariance: We can visualize the spatial 
## and spatio-temporal covariance among species in encounter probability and 
## positive catch rates (depending upon what is turned on via `FieldConfig`):

Cov_List = summarize_covariance(
  Report=Report, 
  ParHat=Obj$env$parList(), 
  Data=TmbData, 
  SD=Opt$SD, 
  plot_cor=FALSE, 
  category_names=levels(Data_Geostat[,'spp']), 
  plotdir=DateFile, 
  plotTF=NULL, 
  mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

## Density surface for each year: We can visualize many types of output from 
## the model. Here I only show predicted density, but other options are 
## obtained via other integers passed to `plot_set` as described in`?plot_maps`

plot_settings = data.frame(plot_num = c(11,1:2,6,7), 
var_name = c('covariates','PoC','PosCPUE', 'epsilon', 'epsilon'), 
stringsAsFactors = F)

for(i in 1:nrow(plot_settings)){

if(!dir.exists(paste0(DateFile, plot_settings$var_name[i] ))) 
  dir.create(paste0(DateFile, plot_settings$var_name[i], '/' ))

plot_maps(plot_set=plot_settings$plot_num[i], 
          Report=Report, 
          TmbData = TmbData,
          Sdreport=Opt$SD, 
          PlotDF=MapDetails_List[["PlotDF"]], 
          MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
          working_dir = paste0(DateFile, plot_settings$var_name[i], '/' ), 
          Year_Set=Year_Set, Years2Include=Years2Include, 
          col = rev(heat.colors(100)),
          mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8,
          category_names=levels(Data_Geostat[,'spp']))
}

if(!dir.exists(paste0(DateFile, 'log_density/' ))) 
  dir.create(paste0(DateFile, 'log_density/' ))

## Plot log-density
plot_maps_density(plot_set=3, 
                  Report=Report, 
                  TmbData = TmbData,
                  Sdreport=Opt$SD, 
                  PlotDF=MapDetails_List[["PlotDF"]], 
                  MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
                  working_dir = paste0(DateFile, 'log_density/' ), 
                  Year_Set=Year_Set, Years2Include=Years2Include, 
                  col = c('white', rev(heat.colors(8))[c(1,4,6,8)]),
                  mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8,
                  category_names=levels(Data_Geostat[,'spp']))

## Index of abundance: The index of abundance is generally most useful for stock 
## assessment models.

Index = plot_biomass_index( DirName=DateFile, 
                            TmbData=TmbData, 
                            Sdreport=Opt[["SD"]], 
                            Year_Set=Year_Set, Years2Include=Years2Include,
                            strata_names=strata.limits[,1], 
                            use_biascorr=TRUE, 
                            category_names=levels(Data_Geostat[,'spp']) )
pander::pandoc.table( Index$Table[,c("Category","Year",
                                     "Estimate_metric_tons","SD_mt")] ) 

## Plot factors: Finally, we can inspect the factor-decomposition for 
## community-level patterns.  This generates many plots, only some of which are 
## included in this tutorial document.

if(!dir.exists(paste0(DateFile, "Factor/"))) dir.create(paste0(DateFile, "Factor/"))

plot_factors( Report=Report, 
              ParHat=Obj$env$parList(), 
              Data=TmbData, 
              SD=Opt$SD, 
              mapdetails_list=MapDetails_List, 
              Year_Set=Year_Set,
              category_names=levels(Data_Geostat[,'spp']), 
              plotdir= paste0(DateFile, "Factor/") )

save(file = paste0(DateFile,"Diagnostics.RData"), 
	list = c('Q', 'Cov_List', 'Enc_prob', 'Index', 'MapDetails_List',
	         'modelno', 'plot_factors', 'plot_residuals', 'plot_settings',
		   'Year_Set', 'Years2Include') )
