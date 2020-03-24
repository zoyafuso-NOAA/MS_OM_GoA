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

figure_wd = c('',
              'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/figure_plot/')[which_machine]

tiff(paste0(figure_wd, 'density_spatiotemporal.tiff'),
     width = 190, height = 190, units = 'mm', res = 200, compression = 'lzw')

par(mar = c(0,0,0,0), mfrow = c(1,3))
for(ispp in 1:3){
  
  
  yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
  plot(1, type = 'n', axes = F, ann = F,
       xlim = range(Extrapolation_List$Data_Extrap[,c('E_km')]),
       ylim = c(min(Extrapolation_List$Data_Extrap[,c('N_km')])-5.55*yrange,
                max(Extrapolation_List$Data_Extrap[,c('N_km')]))
  )
  
  offset = 0
  for(iyear in Years2Include){
    data = as.data.frame(Save$Report$D_gcy[,ispp,iyear])
    names(data) = 'density'
    
    goa = SpatialPointsDataFrame(coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
                                 data = data )
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = 'density')
    
    vals  = Save$Report$D_gcy[,ispp,]
    
    val_cuts = c(1, quantile(vals[vals>1] , na.rm = T, 
                             probs = seq(0,1, length = 9)))
    values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
    
    #offset
    goa_ras = raster::shift(x = goa_ras, dy = -yrange/1.85*offset)
    offset = offset + 1
    
    image(goa_ras, asp = 1, axes = F, ann = F, add = T)
    
    text(x = goa_ras@extent[1] + 0.7*diff(goa_ras@extent[1:2]),
         y = goa_ras@extent[3]+ 0.7*diff(goa_ras@extent[3:4]),
         Year_Set[iyear])
    
  }
  
  mtext(side = 3, Save$Spp[ispp], line = -1.4)
  
  val_cuts = round(val_cuts)
  legend('bottom', fill = hcl.colors(9, "YlOrRd", rev = TRUE), bty = 'n', 
         ncol = 3, cex = 0.90,
         legend = c('<1',paste0(val_cuts[2:(length(val_cuts)-1)], '-',
                                val_cuts[3:length(val_cuts)])) )
  
}

dev.off()
