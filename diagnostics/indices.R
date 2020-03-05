GOA_biomass_indices_wnames = readRDS('C:/Users/zack.oyafuso/Downloads/GOA_biomass_indices_wnames.rds')
spp_df = read.csv("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/spp_df.csv",
                  check.names = F, row.names = 'modelno')

for(spp in colnames(spp_df) ){
  df = subset(GOA_biomass_indices_wnames, 
              SPECIES_NAME == spp & YEAR >= 1996)
  df = df[order(df$YEAR),]
  lower = (df$TOTAL_BIOMASS - sqrt(df$BIOMASS_VAR))/1e6
  upper = (df$TOTAL_BIOMASS + sqrt(df$BIOMASS_VAR))/1e6
  
  plot((TOTAL_BIOMASS/1e6) ~ YEAR, data = df, type = 'b',
       main = spp,
       ylim = c(0, max(upper)), las = 1, 
       ylab = 'Total Biomass (millions of units)')
  
  segments(x0 = df$YEAR, x1 = df$YEAR, 
           y0 = lower, 
           y1 = upper )
  
}

sort(unique(GOA_biomass_indices_wnames$SPECIES_NAME))