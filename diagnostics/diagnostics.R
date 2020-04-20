##################################
## Diagnostic plots
##################################

library(VAST)
library(FishStatsUtils)

modelno = '6j'
factorno = 2

VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output", 
                  modelno, "/Factor_", factorno, "/")
diag_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/diagnostics/VAST_model", 
                  modelno, "/Factor_", factorno, "/")
fun_dir = paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/')
if(! dir.exists(diag_dir) ) dir.create(diag_dir)

#################################
## Load Data and modified functions
#################################
load(paste0(VAST_dir, '/fit.RData'))

for(ifile in c("plot_factors.R", "plot_maps_density.R", "plot_residuals.R",
               "plot_variable_density.R", "summarize_covariance.R")){
  source( paste0(fun_dir, ifile) )
}

load( paste0(dirname(VAST_dir), '/Spatial_Settings_CrVa.RData') )

#We first apply a set of standard model diagnostics to confirm that the model 
## is reasonable and deserves further attention.  If any of these do not look 
## reasonable, the model output should not be interpreted or used.

######################
## Plot data
## It is always good practice to conduct exploratory analysis of data.  Here, I 
## visualize the spatial distribution of data.  Spatio-temporal models involve 
## the assumption that the probability of sampling a given location is 
## statistically independent of the probability distribution for the response 
## at that location.  So if sampling "follows" changes in density, then the 
## model is probably not appropriate!
############################################
Data_Geostat = fit$data_frame
names(Data_Geostat)[c(1:2,5:7)] = c('Lat', 'Lon', 'Catch_KG','Year', "spp")

Spatial_List = fit$spatial_list
Extrapolation_List = fit$extrapolation_list
Opt = fit$parameter_estimates
Obj = fit$tmb_list$Obj
Report = fit$Report
TmbData = fit$data_list


plot_data(Extrapolation_List=Extrapolation_List, 
          Spatial_List=Spatial_List, 
          Data_Geostat=Data_Geostat, 
          PlotDir=diag_dir)

############################################
## Convergence
## Here I print the diagnostics generated during parameter estimation, and I 
## confirm that (1) no parameter is hitting an upper or lower bound and (2) the
## final gradient for each fixed-effect is close to zero. For explanation of 
## parameters, please see `?make_data`.
############################################

pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE',
                                         'Upper','final_gradient')] ) 

max(abs(Opt$diagnostics$final_gradient))

############################################
## Diagnostics for encounter-probability component
## Next, we check whether observed encounter frequencies for either low or high 
## probability samples are within the 95% predictive interval for predicted 
## encounter probability
############################################
Enc_prob = plot_encounter_diagnostic( Report=Report, 
                                      Data_Geostat=Data_Geostat, 
                                      DirName=diag_dir)

############################################
## Diagnostics for positive-catch-rate component
## We can visualize fit to residuals of catch-rates given encounters using a
## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
############################################
if(!dir.exists(paste0(diag_dir,'QQ_Fn/'))) dir.create(paste0(diag_dir,'QQ_Fn/'))
   
Q = plot_quantile_diagnostic( TmbData=TmbData, 
                              Report=Report, 
                              FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", 
                              FileName_Qhist="Q-Q_hist", 
                              DateFile=paste0(diag_dir, 'Q/')) 
 
############################################
## Diagnostics for plotting residuals on a map
## Get region-specific settings for plots
## Decide which years to plot
############################################
MapDetails_List = make_map_info( "Region"='Gulf_of_Alaska', 
                                 "spatial_list"=Spatial_List, 
                                 "Extrapolation_List"=Extrapolation_List )
                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

############################################
## We then plot Pearson residuals.  If there are visible patterns (areas with 
## consistently positive or negative residuals accross or within years) then 
## this is an indication of the model "overshrinking" results towards the 
## intercept, and model results should then be treated with caution.  
############################################
if(!dir.exists(paste0(diag_dir,'Pearson_Residuals/'))) 
  dir.create(paste0(diag_dir,'Pearson_Residuals/'))

plot_residuals(Lat_i=Data_Geostat[,'Lat'], 
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
## Direction of "geometric anisotropy"
## We can visualize which direction has faster or slower decorrelation 
## (termed "geometric anisotropy")
############################################
plot_anisotropy( FileName=paste0(diag_dir,"Aniso.png"), 
                 Report=Report, 
                 TmbData=TmbData )

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
## Density surface for each year
## We can visualize many types of output from the model.  Here I only show 
## predicted density, but other options are obtained via other integers passed 
## to `plot_set` as described in `?plot_maps`
############################################
plot_settings = data.frame(plotno = c(11),
                           foldername = c('covariates'))
irow = 1

if(!dir.exists(paste0(diag_dir, plot_settings$foldername, '/'))) 
  dir.create(paste0(diag_dir, plot_settings$foldername, '/'))

plot_maps(plot_set=plot_settings$plotno[irow],
          Obj = Obj,
          MappingDetails=MapDetails_List[["MappingDetails"]], 
          Report=Report, 
          TmbData = TmbData,
          Sdreport=Opt$SD, 
          PlotDF=MapDetails_List[["PlotDF"]], 
          MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
          Xlim=MapDetails_List[["Xlim"]], 
          Ylim=MapDetails_List[["Ylim"]], 
          FileName=paste0(diag_dir, plot_settings$foldername[irow], '/'), 
          Year_Set=Year_Set, 
          Years2Include=Years2Include, 
          Rotate=MapDetails_List[["Rotate"]],
          Cex=MapDetails_List[["Cex"]], 
          Legend=MapDetails_List[["Legend"]], 
          zone=MapDetails_List[["Zone"]], 
          mar=c(0,0,2,0), 
          oma=c(3.5,3.5,0,0), 
          cex=1.8, 
          category_names=levels(as.factor(Data_Geostat[,'spp'])) )

plot_variable(Y_gt = TmbData$X_gtp[,1],
              map_list = MapDetails_List,
              file_name = 'covariate')


############################################
## Index of abundance
############################################
if(!dir.exists(paste0(diag_dir,'Index/'))) 
  dir.create(paste0(diag_dir,'Index/'))
Index = plot_biomass_index( DirName=paste0(diag_dir,'Index/'), 
                            TmbData=TmbData, 
                            Sdreport=Opt[["SD"]], 
                            Year_Set=Year_Set, 
                            Years2Include=Years2Include, 
                            strata_names=strata.limits[,1], 
                            use_biascorr=TRUE, 
                            category_names=levels(Data_Geostat[,'spp']) )
pander::pandoc.table( Index$Table[,c("Category","Year","Estimate_metric_tons","SD_mt")] ) 

############################################
## Center of gravity and range expansion/contraction
## We can detect shifts in distribution or range expansion/contraction. 
############################################
# plot_range_index(Report=Report, 
#                  TmbData=TmbData, 
#                  Sdreport=Opt[["SD"]], 
#                  Znames=colnames(TmbData$Z_xm), 
#                  PlotDir=diag_dir, 
#                  category_names=levels(Data_Geostat[,'spp']), 
#                  Year_Set=Year_Set)

############################################
## Plot overdispersion
############################################
##We can also plot and inspect overdispersion (e.g., vessel effects, or tow-level fisher targetting), although this example doesn't include any.  

# Plot_Overdispersion( filename1=paste0(diag_dir,"Overdispersion"), 
#                      filename2=paste0(diag_dir,"Overdispersion--panel"), 
#                      Data=TmbData, 
#                      ParHat=ParHat, 
#                      Report=Report, 
#                      ControlList1=list("Width"=5, 
#                                        "Height"=10, 
#                                        "Res"=200, 
#                                        "Units"='in'), 
#                      ControlList2=list("Width"=TmbData$n_c, 
#                                        "Height"=TmbData$n_c, 
#                                        "Res"=200, 
#                                        "Units"='in') )


############################################
## Plot factors
## Finally, we can inspect the factor-decomposition for community-level 
## patterns.  This generates many plots, only some of which are included in 
## this tutorial document.
############################################
if(!dir.exists(paste0(diag_dir,'Factors/'))) 
  dir.create(paste0(diag_dir,'Factors/'))

Plot_factors( Report=Report, 
              ParHat=Obj$env$parList(), 
              Data=TmbData, 
              SD=Opt$SD, 
              mapdetails_list=MapDetails_List, 
              Year_Set=Year_Set, 
              category_names=levels(Data_Geostat$spp), 
              plotdir=paste0(diag_dir,'Factors/') )

############################################
## Save output
############################################
save(list = c('modelno', 'factorno', 'VAST_dir', 'diag_dir', 'fun_dir', 
              'Data_Geostat', 'Spatial_List', 'Extrapolation_List', 'Opt', 
              'Obj', 'Report', 'TmbData', 'Enc_prob', 'Q', 'MapDetails_List', 
              'Year_Set', 'Years2Include', #'Cov_List', 
              'plot_settings', 
              'Index'),
     file = paste0(diag_dir, 'diagnostics.RData'))

