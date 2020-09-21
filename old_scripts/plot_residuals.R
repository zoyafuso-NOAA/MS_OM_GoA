##################################
## Plot Data residuals
##################################

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

modelno = '10a'
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]


VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output", 
                  modelno, "/")

#################################
## Load Data and modified functions
#################################
load(paste0(VAST_dir, '/fit.RData'))

##################################
## Extract objects from fitted object
##################################
Data_Geostat = fit$data_frame
names(Data_Geostat)[c(1:2,5:7)] = c('Lat', 'Lon', 'Catch_KG','Year', "spp")
Data_Geostat$spp = Data_Geostat$spp + 1
Data_Geostat$obs_cpue = Data_Geostat$Catch_KG / Data_Geostat$a_i

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
              "Microstomus pacificus", 
              "Sebastes alutus", "Sebastes B_R", "Sebastes brevispinis", 
              "Sebastes polyspinis", "Sebastes variabilis", 
              "Sebastolobus alascanus" )
ns = length(sci_names)

xrange = range(Extrapolation_List$Data_Extrap[,'E_km'])
yrange = range(Extrapolation_List$Data_Extrap[,'N_km'])
xrange_diff = diff(xrange)
yrange_diff = diff(yrange)

#Grid locations
loc_g = fit$spatial_list$loc_g
loc_i = fit$spatial_list$loc_i

######################################
for(ispp in 1:ns){
  par(mfrow = c(4,3), mar = c(0,0,0,0))
  for(iyear in Years2Include){
    
    #Extract observed data (location and CPUE) for a particular year and spp
    idx = with(Data_Geostat, which(spp == ispp & Year == Year_Set[iyear]) )
    temp_loc_i = loc_i[idx,]
    obs_dens = Data_Geostat$obs_cpue[idx]
    
    #calculate which grids are closest to each withheld data location
    #and extract predicted density
    grid_idx = RANN::nn2(query = temp_loc_i, 
                         data = loc_g, k = 1)$nn.idx
    
    #Extract density of the grids closets
    pred_dens = Report$D_gcy[grid_idx,ispp,iyear]
    
    residual = (obs_dens - pred_dens) / mean(pred_dens)
    abs_residual = log10(abs(residual))
    residual_color = ifelse(residual > 0, 'red', 'blue')
    residual_color = ifelse(abs_residual < 1, 'NA', residual_color)
    residual_size = ifelse(abs_residual > 1, abs_residual, 0)
    
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
      data = data.frame(var = Extrapolation_depths$depth) 
    )
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = 'var')
    
    image(goa_ras, asp = 1, col = rev(terrain.colors(1000)), axes = F, ann = F)
    points(temp_loc_i, pch = 16, 
           cex = residual_size, col = residual_color )
    box()
  } 
  plot(1, type = 'n', axes = F, ann = F)
  text(1,1,sci_names[ispp], cex = 2, font = 3)
}
