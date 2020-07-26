######################################
## Compare Log-Density Estimates
#####################################
rm(list = ls())

#################################
## Import Libraries
#################################
library(VAST)
library(raster)
library(sp)
library(RColorBrewer)
library(plotrix)

#################################
## Set up directories
#################################
modelno_main = '10'

model_log_density = list()
modelno = c('10a','10c')
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]


for(i in 1:2){
  # PP_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/powerpoint_plot/")
  VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output",
                    modelno[i], "/")
  diag_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/diagnostics/",
                    "VAST_model", modelno[i], "/")
  # fun_dir=paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/')
  
  # if(! dir.exists(diag_dir) ) dir.create(diag_dir)
  
  #################################
  ## Load Data and modified functions
  #################################
  load(paste0(VAST_dir, '/fit.RData'))
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
  
  sci_names = c("Atheresthes stomias", "Gadus chalcogrammus", 
                "Gadus macrocephalus", "Glyptocephalus zachirus" , 
                "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
                "Lepidopsetta bilineata", "Lepidopsetta polyxystra",
                "Sebastes brevispinis", "Microstomus pacificus", 
                "Sebastes alutus", "Sebastes B_R", "Sebastes polyspinis", 
                "Sebastes variabilis", "Sebastolobus alascanus" )
  ns = length(sci_names)
  
  xrange = range(Extrapolation_List$Data_Extrap[,'E_km'])
  yrange = range(Extrapolation_List$Data_Extrap[,'N_km'])
  xrange_diff = diff(xrange)
  yrange_diff = diff(yrange)
  
  model_log_density[[i]] = Report$D_gcy
}

{
  for(ispp in 1:ns){
    # png(paste0(diag_dir, 'Density/density_', 
    #            gsub(sci_names[ispp], pattern = ' ', replacement = '_'), 
    #            '.png'), 
    #     units = 'in', height = 5, width = 7, res = 500)
    par(mar = c(0,0,0,0), oma = rep(0.5,4))
    par(mfrow = c(4,3))
    
    for(iyear in Years2Include){
      vals = log(model_log_density[[2]][,ispp,iyear]/model_log_density[[1]][,ispp,iyear])
      
      val_cuts = log(c(1/(10:1),  2:10))
      
      goa = SpatialPointsDataFrame(
        coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
        data = data.frame(density = vals) )
      goa_ras = raster(goa, resolution = 5)
      goa_ras = rasterize(x = goa, y = goa_ras, field = 'density')
      
      values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
      image(goa_ras,
            asp = 1, axes = F, ann = F, 
            col = colorRampPalette(c('red', 'white', 'darkblue'))(length(val_cuts)-1))
      
      text(x = goa_ras@extent[1] + 0.7*diff(goa_ras@extent[1:2]),
           y = goa_ras@extent[3]+ 0.7*diff(goa_ras@extent[3:4]),
           Year_Set[iyear], cex = 1)
      
      # val_cuts = round(val_cuts[-1])
      
      # legend('bottom', fill = colors, bty = 'n',
      #        ncol = 3, cex = 0.65,
      #        legend = c('<1', paste0('1-', val_cuts[2]), 
      #                   paste0(val_cuts[2:(length(val_cuts)-1)], '-',
      #                          val_cuts[3:length(val_cuts)])) )
      box()  
    }
    plot(1, type = 'n', axes = F, ann = F)
    text(1,1, sci_names[ispp], cex = 1.5, font = 3)
    # plotrix::color.legend(xl = 0.65,
    #                       xr = 1.36,
    #                       yb = 0.75,
    #                       yt = 0.85,
    #                       legend = pretty(c(-max_val, max_val)),
    #                       rect.col = colorRampPalette(c('red',
    #                                                     'white',
    #                                                     'darkblue'))(9),
    #                       gradient = 'x', align = 'rb')
    
    # dev.off()
  }  
}
