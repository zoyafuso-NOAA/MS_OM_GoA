###############################################################################
## Project:       Compare Design- and VAST-based abundance indices
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot mean and CIs of abundance indices
###############################################################################

##################################################
####   Set up directories
##################################################
rm(list = ls())

modelno = '10b'
VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/", 
                  "VAST_output", modelno, "/")
diag_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/diagnostics/",
                  "VAST_model", modelno, "/")
github_dir = paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/')

##################################################
#### Import: Species Set
#### Design-based estimates
#### VAST fit
##################################################
spp_df = read.csv(paste0(github_dir, "data/spp_df.csv"), 
                  check.names=F, header = T, row.names = 'modelno')

GOA_DBE = readRDS(file = paste0("C:/Users/Zack Oyafuso/Documents/GitHub/",
                                "MS_OM_GoA/data/",
                                "GOA_biomass_indices_wnames.rds") )

load(file = paste0(VAST_dir, 'fit.RData'))

##################################################
#### Constants
##################################################
Year_Set = seq(min(fit$data_frame[,'t_i']),max(fit$data_frame[,'t_i']))
Years2Include = which( Year_Set %in% sort(unique(fit$data_frame[,'t_i'])))

Index_Ests = fit$Report$Index_cyl[,Years2Include,1]
Index_SDs = matrix(data = fit$parameter_estimates$SD$sd[attributes(fit$parameter_estimates$SD$value)$names == "Index_cyl"], nrow = 15 )[,Years2Include]

sci_names = names(spp_df)[unlist(spp_df[modelno,])] 
sci_names = sort(sci_names)
spp_to_include = which(!(sci_names %in% 'Sebastes B_R') == T)
ns = length(spp_to_include)

{
  png(paste0(diag_dir, 'indices_comparison_DBE.png'), width = 12, height = 6, 
      units = 'in', res = 500)
  par(mfrow = c(3,5), mar = c(3,3,2,1), oma = c(0,2.5,0,0))
  for(spp in spp_to_include) {
    
    #Design based Estimator and SD Interval
    temp_DBE = subset(GOA_DBE, SPECIES_NAME == sci_names[spp] &
                        YEAR %in% Year_Set[Years2Include])
    
    temp_DBE = temp_DBE[order(temp_DBE$YEAR),]
    
    upper_DBE = (temp_DBE$TOTAL_BIOMASS + sqrt(temp_DBE$BIOMASS_VAR))/1e6
    lower_DBE = (temp_DBE$TOTAL_BIOMASS - sqrt(temp_DBE$BIOMASS_VAR))/1e6
    
    #VAST SD Intervals
    upper = (Index_Ests + Index_SDs)/1e6
    lower = (Index_Ests - Index_SDs)/1e6
    
    #Empty plot
    plot(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6, type = 'n', 
         ylab = 'Index', xlab = "Year", las = 1, 
         ylim = c(0, max(c(upper[spp,], upper_DBE) )) )
    
    #Plot VAST interals
    polygon(x = c(Year_Set[Years2Include],
                  rev(Year_Set[Years2Include])),
            y = c(lower[spp,], rev(upper[spp,])), col = 'grey', lty = 'dotted')
    
    lines(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6)
    points(x = Year_Set[Years2Include], Index_Ests[spp,]/1e6, pch= 16)
    
    #Plot design-based intervals
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
  
  #Plot Legend
  plot(1, type = 'n', axes = F, ann = F)
  legend('center', legend = c('DBE', 'VAST'), 
         col = c('red', 'black'), pch = 16, lty=1, cex = 3)
  mtext(side = 2, outer = T, text = 'Abundance Index (million metric tons)',
        line = 1)
  dev.off()
}

