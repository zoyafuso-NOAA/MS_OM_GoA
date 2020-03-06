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

frame1 <- buildFrameDF(df = df,
                       id = "x",
                       X = c("depth"),
                       Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                       domainvalue = "Domain")
strata1 <- buildStrataDF(frame1, progress=F)




# checkInput(errors = checkInput(errors = cv, 
#                                strata = strata1, 
#                                sampframe = frame1))
# 
# allocation <- bethel(strata1,cv[1,])
# sum(allocation)

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = "depth",
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

#Settings for optimizer
settings = results = expand.grid(cv = c(0.2),
                                 pops = c(25,50,100),
                                 minnumstr = c(25, 50),
                                 mut_change = c(0.01, 0.05, 0.1, 0.5),
                                 elitism_rate = c(0.1, 0.2, 0.5))

ns = Save$TmbData$n_c

# rm(list  = ls()[!ls() %in% c('settings', 'frame', 'modelno', 'ns') ])

library(foreach)
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

foreach(i = 1:nrow(settings), 
        .packages = 'SamplingStrata',
        .export = c('settings', 'frame',  'modelno', 'ns') ) %dopar% 
  {
    library(SamplingStrata)
    wd = paste0("C:/Users/Zack Oyafuso/Documents/",
                "GitHub/MS_OM_GoA/Optimum_Allocation/",
                "model_", modelno, "/",
                'cv_', settings$cv[i], '_', 
                'pop_', settings$pops[i], '_',
                'minnumstr_', settings$minnumstr[i], '_',
                'mutchange_', settings$mut_change[i], '_',
                'elitism_rate_', settings$elitism_rate[i], '.RData')
    
    cv = list()
    for(i in 1:ns) cv[[paste0('CV', i)]] = settings$cv[i]
    cv[['DOM']] = 'GoA'; cv[['domainvalue']]=1
    cv <- as.data.frame(cv)
    # cv
    
    init_sol <- KmeansSolution2(frame=frame,
                                errors=cv,
                                maxclusters = 10)  
    
    nstrata <- tapply(init_sol$suggestions,
                      init_sol$domainvalue,
                      FUN=function(x) length(unique(x)))
    # nstrata
    
    initial_solution <- prepareSuggestion(init_sol,frame,nstrata)
    
    set.seed(1234 + i)
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 50,
                            pops = settings$pops[i],
                            minnumstr=settings$minnumstr[i],
                            mut_chance = settings$mut_change[i],
                            nStrata = nstrata,
                            suggestions = initial_solution,
                            showPlot = F,
                            parallel = F)
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress=FALSE)
    
    save(list=c('strataStructure', 'solution'), file = wd)
  }



