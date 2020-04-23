##################################
## Diagnostic plots
##################################

library(VAST)
library(raster)
# library(FishStatsUtils)

modelno_main = '7'
modelno = '7a'

VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output", 
                  modelno_main, "/VAST_output", modelno, "/")
diag_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/diagnostics/",
                  "VAST_model", modelno_main, "/VAST_output", modelno, "/")
fun_dir = paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/')
if(! dir.exists(diag_dir) ) dir.create(diag_dir)

#################################
## Load Data and modified functions
#################################
load(paste0(VAST_dir, '/fit.RData'))

for(ifile in c("plot_factors.R", "plot_maps_density.R", "plot_residuals.R",
               "plot_variable_density.R", "summarize_covariance.R")){
  source( paste0(fun_dir, ifile) )
  rm(ifile)
}

load( paste0(dirname(VAST_dir), '/Spatial_Settings_CrVa.RData') )

##################################
## Extract objects from fitted object
##################################
Data_Geostat = fit$data_frame
names(Data_Geostat)[c(1:2,5:7)] = c('Lat', 'Lon', 'Catch_KG','Year', "spp")

Spatial_List = fit$spatial_list
Extrapolation_List = fit$extrapolation_list
Opt = fit$parameter_estimates
Obj = fit$tmb_list$Obj
Report = fit$Report
TmbData = fit$data_list

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

MapDetails_List = make_map_info( "Region"='Gulf_of_Alaska', 
                                 "spatial_list"=Spatial_List, 
                                 "Extrapolation_List"=Extrapolation_List )

sci_names = levels(
  read.csv(paste0(dirname(fun_dir), 
                  '/data/data/GOA_multspp.csv'))$SPECIES_NAME) 
ns = length(sci_names)

############################################
## Diagnostics for positive-catch-rate component
## We can visualize fit to residuals of catch-rates given encounters using a
## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
############################################
   
Q = plot_quantile_diagnostic( TmbData=TmbData, 
                              Report=Report, 
                              FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", 
                              FileName_Qhist="Q-Q_hist", 
                              DateFile=diag_dir ) 

############################################
## We then plot Pearson residuals.  If there are visible patterns (areas with 
## consistently positive or negative residuals accross or within years) then 
## this is an indication of the model "overshrinking" results towards the 
## intercept, and model results should then be treated with caution.  
############################################
if(!dir.exists(paste0(diag_dir,'Pearson_Residuals/'))) 
  dir.create(paste0(diag_dir,'Pearson_Residuals/'))

PResid = plot_residuals(Lat_i=Data_Geostat[,'Lat'], 
                        Lon_i=Data_Geostat[,'Lon'], 
                        TmbData=TmbData, 
                        Report=Report, 
                        Q=Q, 
                        spatial_list = Spatial_List,
                        extrapolation_list = Extrapolation_List,
                        working_dir=paste0(diag_dir,'Pearson_Residuals/'), 
                        Year_Set=Year_Set, 
                        Years2Include=Years2Include, 
                        mar=c(0,0,2,0), 
                        oma=c(3.5,3.5,0,0), 
                        cex=1.8)


############################################
## Plot spatial and spatio-temporal covariance
## We can visualize the spatial and spatio-temporal covariance among species 
## in encounter probability and positive catch rates (depending upon what is
### turned on via `FieldConfig`):
############################################
Cov_List = summarize_covariance( Report=Report, 
                                 ParHat=Obj$env$parList(), 
                                 Data=TmbData, 
                                 SD=Opt$SD, 
                                 plot_cor=FALSE, 
                                 category_names=levels(Data_Geostat[,'spp']), 
                                 plotdir=diag_dir, 
                                 mgp=c(2,0.5,0), 
                                 tck=-0.02, 
                                 oma=c(0,10,2,2) )

############################################
## Save output
############################################
save(list = c('Data_Geostat', 'Spatial_List', 'Extrapolation_List', 'Opt', 
              'Obj', 'Report', 'TmbData', 'Q', 'MapDetails_List', 'Year_Set', 
              'Years2Include', 'sci_names', 'ns', 'PResid', 'Cov_List'),
     file = paste0(diag_dir, 'diagnostics.RData'))

