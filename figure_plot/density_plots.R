################################
## General Plotting Function
################################

rm(list = ls())
which_machine = c('Zack_MAC' = 1)
modelno = '6g'
results_dir = paste0(c('/Users/zackoyafuso/Google Drive/')[which_machine],
                     'VAST_Runs/VAST_Output',  modelno)
setwd(results_dir)

load('VAST_MS_GoA_Run.RData')
load('Spatial_Settings.RData')

MapDetails_List = make_map_info( "Region"='Gulf_of_Alaska', 
                                 "spatial_list"=Spatial_List, 
                                 "Extrapolation_List"=Extrapolation_List )
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

for(iyear in Years2Include){
 goa = SpatialPointsDataFrame(coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
                              data = data.frame(X1=res_df[,winner]) )
}

