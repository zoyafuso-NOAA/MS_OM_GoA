rm(list = ls())

GOA_DBE = readRDS(file = "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/GOA_biomass_indices_wnames.rds")

setwd("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/")
modelno = '6c'

load(file = paste0('VAST_output', modelno, '/', 'VAST_MS_GoA_Run.RData'))
load(file = paste0('VAST_output', modelno, '/', 'Spatial_Settings.RData'))

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Index_Ests = Save$Report$Index_cyl[,Years2Include,1]
Index_SDs = matrix(data = Save$Opt$SD$sd[attributes(Save$Opt$SD$value)$names == "Index_cyl"], nrow = length(Save$Spp))[,Years2Include]


{tiff(paste0('diagnostics/VAST_model', modelno, '/indices_comparison_DBE.tiff'),
     width = 6, height = 6, units = 'in', res = 500, compression = 'lzw')
par(mfrow = c(5,3), mar = c(1,3,2,1))
for(spp in (1:length(Save$Spp))[-12] ){
  
  temp_DBE = subset(GOA_DBE, SPECIES_NAME == Save$Spp[spp] & YEAR %in% Year_Set[Years2Include])
  
  temp_DBE = temp_DBE[order(temp_DBE$YEAR),]
  
  upper_DBE = (temp_DBE$TOTAL_BIOMASS + sqrt(temp_DBE$BIOMASS_VAR))/1e6
  lower_DBE = (temp_DBE$TOTAL_BIOMASS - sqrt(temp_DBE$BIOMASS_VAR))/1e6
  
  upper = (Index_Ests + Index_SDs)/1e6
  lower = (Index_Ests - Index_SDs)/1e6
  
  plot(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6, type = 'n', ylab = 'Index', xlab = "Year",las = 1, ylim = c(0, max(c(upper[spp,], upper_DBE) )),
       main = )
  
  polygon(x = c(Year_Set[Years2Include],
                rev(Year_Set[Years2Include])),
          y = c(lower[spp,], rev(upper[spp,])), col = 'grey', lty = 'dotted')
  
  lines(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6)
  points(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6, pch= 16)
  
  points(Year_Set[Years2Include],
         temp_DBE$TOTAL_BIOMASS/1e6, col = 'red', pch = 16)
  lines(Year_Set[Years2Include],
        temp_DBE$TOTAL_BIOMASS/1e6, col = 'red')
  
  segments(x0 = Year_Set[Years2Include],
           x1 = Year_Set[Years2Include],
           y0 = lower_DBE,
           y1 = upper_DBE, 
           col = 'red')
  
  mtext(side = 3, unique(temp_DBE$COMMON_NAME), font = 1)
  
}

dev.off()
}
