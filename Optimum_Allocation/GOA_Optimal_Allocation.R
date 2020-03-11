#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata); library(sp)

VAST_model = "6c"
setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

wd = paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
            "Optimum_Allocation/model_", VAST_model)

if(!dir.exists(wd)) dir.create(wd)

load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))
load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

df = cbind(
  data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                          breaks = c(-171, -159, -154, -147, -140, -130), 
                          labels = c('Shumagin_1', 'Chirikof_2', 'Kodiak_3',
                                     'Yakutak_4', 'SE_5')),
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X = Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean)
)
names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth"),#, 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.3),
                       mut_change = c(0.01, 0.1, 0.5),
                       elitism_rate = c(0.1, 0.2, 0.5))

ns = Save$TmbData$n_c
domains = unique(df$Domain)
ndom = length(unique(frame$domainvalue))
  
# rm(list  = ls()[!ls() %in% c('settings', 'frame', 'VAST_model', 'ns') ])

# for(i in 1:nrow(settings)){
#   
#   wd = paste0("C:/Users/Zack Oyafuso/Documents/",
#               "GitHub/MS_OM_GoA/Optimum_Allocation/",
#               "model_", VAST_model, "/",
#               'cv_', settings$cv[i], '_', 
#               'mut_change_', settings$mut_change[i], '_',
#               'elitism_rate_', settings$elitism_rate[i], '.RData')
  i=1
  cv = list()
  for(spp in 1:ns) cv[[paste0('CV', spp)]] = rep(settings$cv[i], ndom)
  cv[['DOM']] = levels(domains)
  cv[['domainvalue']] = as.numeric(domains)
  cv <- as.data.frame(cv)
  
  set.seed(1234 + i)
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 20,
                          pops = 50,
                          elitism_rate = settings$elitism_rate[i],
                          mut_chance = settings$mut_change[i],
                          nStrata = rep(3, ndom),
                          showPlot = T,
                          parallel = T)
  
  strataStructure <- summaryStrata(solution$framenew,
                                   solution$aggr_strata,
                                   progress=FALSE)
  
#   save(list=c('strataStructure', 'solution'), file = wd)
# }

  library(RColorBrewer); library(raster)
  load('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Extrapolation_depths.RData')
  
  goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                               data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                            X1 = solution$indices$X1,
                                            domain = df$Domain) )
  # goa = SpatialPoints(coords = Extrapolation_depths[,c('E_km', 'N_km')] )
  goa_ras = raster(goa, resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'X1')
  
  par(mar = c(0,0,0,0))
  nstrata = nrow(strataStructure)
  plot(goa_ras, col = rev(brewer.pal(n = 3, name= 'Spectral')), legend = T )
  
  
  # plot(N_km ~ E_km, data = Extrapolation_depths, pch = '.', cex = 2,
       # col = brewer.pal(n = 5, name = 'Spectral')[frame$domainvalue])
  
  expected_CV(solution$aggr_strata)
  sum(strataStructure$Allocation)

  