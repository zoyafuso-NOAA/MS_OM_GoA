library(raster); library(RColorBrewer); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')
load('../Extrapolation_depths.RData')

settings = expand.grid(cv = c(0.1, 0.3),
                       mut_change = c(0.01, 0.1, 0.5),
                       elitism_rate = c(0.1, 0.2, 0.5))
VAST_model = '6'

for(i in 1:nrow(settings)){
  wd = paste0("C:/Users/Zack Oyafuso/Documents/",
              "GitHub/MS_OM_GoA/Optimum_Allocation/",
              "model_", VAST_model, "/",
              'cv_', settings$cv[i], '_', 
              'mut_change_', settings$mut_change[i], '_',
              'elitism_rate_', settings$elitism_rate[i], '.RData')
  
  load(wd)
  settings$n[i] = sum( strataStructure$Allocation ) 
}

winner = which.min(settings$n)
wd = paste0("C:/Users/Zack Oyafuso/Documents/",
            "GitHub/MS_OM_GoA/Optimum_Allocation/",
            "model_", VAST_model, "/",
            'cv_', settings$cv[winner], '_', 
            'nstrata_', settings$nstrata[winner], '_',
            'mut_change_', settings$mut_change[winner], '_',
            'elitism_rate_', settings$elitism_rate[winner], '.RData')

load(wd)

ns = sum( substr(x = colnames(solution$framenew), start = 1, stop = 1) == 'Y')
goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                             data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                        X1 = solution$indices$X1) )
# goa = SpatialPoints(coords = Extrapolation_depths[,c('E_km', 'N_km')] )
goa_ras = raster(goa, resolution = 5)
goa_ras =rasterize(x = goa, y = goa_ras, field = 'X1')

par(mar = c(0,0,0,0))
nstrata = nrow(strataStructure)
plot(goa_ras, col = brewer.pal(n = nstrata, name= 'Paired'), legend = T )

expected_CV(solution$aggr_strata)
sum(strataStructure$Allocation)

