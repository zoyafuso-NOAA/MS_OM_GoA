#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata); 
library(sp); library(raster)

VAST_model = "6g"
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

yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))


figure_wd = c('',
              'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/figure_plot/')[which_machine]


tiff(paste0(figure_wd, 'Mean_CV.tiff'),
     width = 190, height = 200, units = 'mm', res = 200, compression = 'lzw')
par(mar = c(0,0,0,0), mfrow = c(1,2))

for(itype in c('Mean', 'CV')){
    plot(1, type = 'n', axes = F, ann = F,
         xlim = range(Extrapolation_List$Data_Extrap[,c('E_km')]),
         ylim = c(min(Extrapolation_List$Data_Extrap[,c('N_km')])-9*yrange,
                  max(Extrapolation_List$Data_Extrap[,c('N_km')]))
    )
    
    offset = 0
    
    
    for(i in 1:length(Save$Spp)){
        
        temp = as.data.frame(cbind(Ests_means[,i], Ests_CV[,i]))
        names(temp) = c("Mean", 'CV')
        
        goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 
                                                                      'N_km')], 
                                     data = temp)
        
        goa_ras = raster(goa, resolution = 5)
        goa_ras =rasterize(x = goa, y = goa_ras, field = itype)
        
        vals  = values(goa_ras)
        
        if(itype == 'CV') val_cuts = c(0,0.25,0.5,0.75,1,2,10)
        
        if(itype == 'Mean')     {
            vals = vals[vals > 1]
            val_cuts = c(0, quantile(vals , na.rm = T, 
                                     probs = seq(0,1, length = 11))[-1])
        } 
        
        values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
        
        #offset
        goa_ras = raster::shift(x = goa_ras, dy = -yrange/1.5*offset)
        offset = offset + 1
        
        image(goa_ras, asp = 1, axes = F, ann = F, add = T,
              col = hcl.colors(c('Mean' = 10, 'CV' = 6)[itype], 
                               "YlOrRd", rev = TRUE))
        
        text(x = goa_ras@extent[1] + 0.175*diff(goa_ras@extent[1:2]),
             y = goa_ras@extent[3]+ 0.35*diff(goa_ras@extent[3:4]),
             Save$Spp[i], cex = 0.6, srt = 10)
        
    }
    mtext(side = 3, line = -1, text = itype)
}
dev.off()
