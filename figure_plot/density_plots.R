################################
## General Plotting Function
################################
rm(list = ls())

library(VAST); library(sp); library(raster)

which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
modelno = '6g'
results_dir = paste0(c('/Users/zackoyafuso/Google Drive/',
                       'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                     'VAST_Runs/VAST_Output',  modelno)
setwd(results_dir)

load('VAST_MS_GoA_Run.RData')
load('Spatial_Settings.RData')

MapDetails_List = make_map_info( "Region"='Gulf_of_Alaska', 
                                 "spatial_list"=Spatial_List, 
                                 "Extrapolation_List"=Extrapolation_List )
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

par(mar = c(0,0,0,0))
yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
plot(1, type = 'n', 
     xlim = range(Extrapolation_List$Data_Extrap[,c('E_km')]),
     ylim = c(min(Extrapolation_List$Data_Extrap[,c('N_km')])-3*yrange,
              max(Extrapolation_List$Data_Extrap[,c('N_km')]))
     )

offset = 0
ispp = 1
for(iyear in Years2Include){
  data = as.data.frame(Save$Report$D_gcy[,ispp,iyear])
  names(data) = 'density'
  
  goa = SpatialPointsDataFrame(coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
                                 data = data )
  goa_ras = raster(goa, resolution = 5)
  goa_ras = rasterize(x = goa, y = goa_ras, field = 'density')
  
  val_cuts = quantile( Save$Report$D_gcy[,ispp,], na.rm = T )
  values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
  
  #offset
  goa_ras = raster::shift(x = goa_ras, dy = -yrange/1.75*offset)
  offset = offset + 1

  image(goa_ras, asp = 1, axes = F, ann = F, add = T)
  
  # val_cuts = round(val_cuts)
  # legend('bottomright', fill = hcl.colors(4, "YlOrRd", rev = TRUE), bty = 'n',
  #        legend = paste0(val_cuts[1:(length(val_cuts)-1)], '-', 
  #                                       val_cuts[2:length(val_cuts)])) 
}

