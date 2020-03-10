#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata); 
library(sp); library(raster)

VAST_model = "6c"
setwd('C:/Users/zack.oyafuso/Desktop/VAST_Runs/')

load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))
load("C:/Users/zack.oyafuso/Work/Github/MS_OM_GoA/Extrapolation_depths.RData")
# load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Ests = Save$Report$Index_gcyl[,,Years2Include,1]
TotVar = apply(X = Ests, MARGIN = c(1,2), FUN = function(x) var(as.vector(x)) )
Ests_means =apply(Ests, MARGIN = 1:2, mean)
Ests_CV = sqrt(TotVar)/ Ests_means

par(mar = c(1,1,4,4), mfrow = c(1,2))
for(i in 1:length(Save$Spp)){
    temp = as.data.frame(cbind(Ests_means[,i], Ests_CV[,i]))
    names(temp) = c("Mean", 'CV')
    
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                                 data = temp)
    
    goa_ras = raster(goa, resolution = 5)
    
    goa_ras2 =rasterize(x = goa, y = goa_ras, field = 'Mean')
    
    plot(cut(goa_ras2, quantile(values(goa_ras2), na.rm = T)), main = paste(Save$Spp[i], 'Mean'), axes = F, col = rev(terrain.colors(4)) )
    
    goa_ras2 =rasterize(x = goa, y = goa_ras, field = 'CV')
    
    plot(cut(goa_ras2, quantile(values(goa_ras2), na.rm = T)), main = paste(Save$Spp[i], 'CV'), axes = F, col = rev(terrain.colors(4)) )
    
}


# goa = SpatialPoints(coords = Extrapolation_depths[,c('E_km', 'N_km')] )

par(mar = c(0,0,2,0), mfrow = c(4,4))
for(spp in Save$Spp){

}

