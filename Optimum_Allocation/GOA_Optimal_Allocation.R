#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################


rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

VAST_model = "2b"
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

cv = list()
for(i in 1:Save$TmbData$n_c) cv[[paste0('CV', i)]] = 0.1
cv[['DOM']] = 'GoA'; cv[['domainvalue']]=1
cv <- as.data.frame(cv)
cv


checkInput(errors = checkInput(errors = cv, 
                               strata = strata1, 
                               sampframe = frame1))

allocation <- bethel(strata1,cv[1,])
sum(allocation)

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = "depth",
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

init_sol <- KmeansSolution2(frame=frame,
                            errors=cv,
                            maxclusters = 10)  

nstrata <- tapply(init_sol$suggestions,
                  init_sol$domainvalue,
                  FUN=function(x) length(unique(x)))
nstrata

#Settings for optimizer
settings = expand.grid(pops = c(25,50,100),
                       minnumstr = c(10, 20),
                       mut_change = c(0.01, 0.05, 0.1),
                       elitism_rate = c(0.1, 0.2, 0.5))

library(foreach)
initial_solution <- prepareSuggestion(init_sol,frame,nstrata)

rm(list  = ls()[!ls() %in% c('settings', 'cv', 'frame', 
                             'nstrata', 'initial_solution', 
                             'modelno') ])

cl <- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

foreach(i = 1:nrow(settings), 
        .packages = 'SamplingStrata',
        .export = c('settings', 'cv', 'frame', 
                    'nstrata', 'initial_solution', 'modelno') ) %do% 
  {
    library(SamplingStrata)
    wd = paste0("C:/Users/Zack Oyafuso/Documents/",
                "GitHub/MS_OM_GoA/Optimum_Allocation/model_", modelno, "/",
                'pop_', settings$pops[i], '_',
                'minnumstr_', settings$minnumstr[i], '_',
                'mutchange_', settings$mut_change[i], '_',
                'elitism_rate_', settings$elitism_rate[i], '.RData')
    
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
                            showPlot = T,
                            parallel = F)
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress=FALSE)
    
    save(list=c('strataStructure', 'solution'), file = wd)
  }



