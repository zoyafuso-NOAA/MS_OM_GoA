library(raster); library(RColorBrewer); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')
load('../Extrapolation_depths.RData')

settings = results = expand.grid(cv = c(0.2),
                                 pops = c(25,50,100),
                                 minnumstr = c(25, 50),
                                 mut_change = c(0.01, 0.05, 0.1, 0.5),
                                 elitism_rate = c(0.1, 0.2, 0.5))
modelno = '4b'

for(i in 1:nrow(settings)){
  wd = paste0("model_", modelno, "/",
              'cv_0.15_', 
              'pop_', settings$pops[i], '_',
              'minnumstr_', settings$minnumstr[i], '_',
              'mutchange_', settings$mut_change[i], '_',
              'elitism_rate_', settings$elitism_rate[i], '.RData')
  
  load(wd)
  settings$n[i] = sum( strataStructure$Allocation ) 
}

winner = which.min(settings$n)
wd = paste0("model_", modelno, "/",
            'cv_0.15_', 
            'pop_', settings$pops[winner], '_',
            'minnumstr_', settings$minnumstr[winner], '_',
            'mutchange_', settings$mut_change[winner], '_',
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

