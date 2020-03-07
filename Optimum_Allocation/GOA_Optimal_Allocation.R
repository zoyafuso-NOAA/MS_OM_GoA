#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

VAST_model = "4b"
load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))
load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

df = cbind(
  data.frame(Domain = "GoA",
             x = 1:Save$TmbData$n_g,
             depth = Extrapolation_depths$depth),
  apply(X = Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean)
)
names(df)[-(1:3)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = "depth",
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.1, 0.2, 0.3),
                       nstrata = c(5,10),
                       mut_change = c(0.01, 0.1, 0.5),
                       elitism_rate = c(0.1, 0.2, 0.5))

ns = Save$TmbData$n_c

rm(list  = ls()[!ls() %in% c('settings', 'frame', 'modelno', 'ns') ])

library(foreach)
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

foreach(i = 1:nrow(settings), 
        .packages = 'SamplingStrata' ) %dopar% 
  {
    library(SamplingStrata)
    wd = paste0("C:/Users/Zack Oyafuso/Documents/",
                "GitHub/MS_OM_GoA/Optimum_Allocation/",
                "model_", modelno, "/",
                'cv_', settings$cv[i], '_', 
                'nstrata_', settings$nstrata[i], '_',
                'mut_change_', settings$mut_change[i], '_',
                'elitism_rate_', settings$elitism_rate[i], '.RData')
    
    cv = list()
    for(i in 1:ns) cv[[paste0('CV', i)]] = settings$cv[i]
    cv[['DOM']] = 'GoA'; cv[['domainvalue']] = 1
    cv <- as.data.frame(cv)
    
    set.seed(1234 + i)
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 50,
                            pops = 100,
                            elitism_rate = settings$elitism_rate[i],
                            mut_chance = settings$mut_change[i],
                            nStrata = settings$nstrata[i],
                            showPlot = F,
                            parallel = F)
    
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress=FALSE)
    
    save(list=c('strataStructure', 'solution'), file = wd)
  }



