#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata); 
library(sp); library(raster)

VAST_model = "6c"
setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')

load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))
# load("C:/Users/zack.oyafuso/Work/Github/MS_OM_GoA/Extrapolation_depths.RData")
load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Ests = Save$Report$Index_gcyl[,,Years2Include,1]
TotVar = apply(X = Ests, MARGIN = c(1,2), FUN = function(x) var(as.vector(x)) )
Ests_means =apply(Ests, MARGIN = 1:2, mean)
Ests_CV = sqrt(TotVar)/ Ests_means




for(i in 1:length(Save$Spp)){
    if(i %in% c(1, 5, 9, 13)){
        plot_number = ceiling(i/4)
        plot_name = paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA",
                           "/Optimum_Allocation/MeanCV_plots", plot_number,
                           ".tiff")
        tiff(filename = plot_name, width = 12, height =6, units = 'in', 
             res = 500, compression = 'lzw')
        par(mar = c(2,2,2,4), mfrow = c(2,4), oma = c(1,0,0,0))
    }
    temp = as.data.frame(cbind(Ests_means[,i], Ests_CV[,i]))
    names(temp) = c("Mean", 'CV')
    
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 
                                                                  'N_km')], 
                                 data = temp)
    
    goa_ras = raster(goa, resolution = 5)
    
    goa_ras2 =rasterize(x = goa, y = goa_ras, field = 'Mean')
    
    goa_ras_plot = cut( goa_ras2, quantile(values(goa_ras2), na.rm = T))
    
    plot(goa_ras_plot, 
         axes = F, col = rev(terrain.colors(4)), legend = F )
    mtext(side = 3, paste(Save$Spp[i], 'Mean'), cex =  0.75)
    
    legend_x = c(1.05, 1.1); legend_y = c(0.05, 0.85)
    xl = (1 - legend_x[1]) * par("usr")[1] + (legend_x[1]) * 
        par("usr")[2]
    xr = (1 - legend_x[2]) * par("usr")[1] + (legend_x[2]) * 
        par("usr")[2]
    yb = (1 - legend_y[1]) * par("usr")[3] + (legend_y[1]) * 
        par("usr")[4]
    yt = (1 - legend_y[2]) * par("usr")[3] + (legend_y[2]) * 
        par("usr")[4]
    align = c("lt", "rb")[2]
    gradient = c("x", "y")[2]
    
    plotrix::color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                          legend = round(quantile(values(goa_ras2), 
                                                  seq(0.25, 1, 0.25), 
                                                  na.rm = T), 2), 
                          rect.col = rev(terrain.colors(4)), 
                          cex = 0.75, 
                          align = align, 
                          gradient = gradient)

goa_ras2 =rasterize(x = goa, y = goa_ras, field = 'CV')

plot(cut(goa_ras2, quantile(values(goa_ras2), na.rm = T)), axes = F, col = rev(terrain.colors(4)), legend = F )

mtext(side = 3, paste(Save$Spp[i], 'CV'), cex = 0.75)

plotrix::color.legend(xl = xl, yb = yb, xr = xr, yt = yt, 
                      legend = round(quantile(values(goa_ras2), 
                                              seq(0.25, 1, 0.25), 
                                              na.rm = T), 2), 
                      rect.col = rev(terrain.colors(4)), 
                      cex = 0.75, 
                      align = align, 
                      gradient = gradient)

if(i %in% c(4,8,12,15)) dev.off()
}


