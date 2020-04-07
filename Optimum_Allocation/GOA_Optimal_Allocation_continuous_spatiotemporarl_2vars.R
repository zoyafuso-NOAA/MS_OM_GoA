################################################
## Optimization: Lowest CV for a given sample size
################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(VAST);  library(mvtnorm); library(sp); library(RColorBrewer); 
library(raster)
library(memoise); library(doParallel); library(foreach); library(iterators); 
library(parallel); library(pbapply); library(formattable)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[4]

SamplingStrata_dir = paste0(c('', 
                              'C:/Users/Zack Oyafuso',
                              'C:/Users/zack.oyafuso',
                              'C:/Users/zack.oyafuso')[which_machine],
                            '/Downloads/SamplingStrata-master/R')
github_dir = paste0(c('', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')
VAST_model = "6g"
VAST_dir = paste0(c('', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model)

output_wd = c(paste0('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
                     'Optimum_Allocation/model_', VAST_model),
              paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model))[which_machine]

#########################
## Load functions from SamplingStrata packages into global environment
## Load modified buildStrataDF function
#########################
for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(github_dir, '/buildStrataDF_Zack.R'))

###########################
## Load Data
###########################
load(paste0(output_wd, '/optimization_data_model_', 
            VAST_model, '.RData'))

############################
## Settings for optimizer
############################
nstrata = unlist(list('Zack_MAC'=NA, 
               'Zack_PC' =c(5,7,10, 15), 
               'Zack_GI_PC'= c(20,25,30), 
               'VM' = c(40,50,60))[which_machine])

nstrata = 40

settings = data.frame()
res_df = data.frame(id = 1:nrow(frame))
strata_list = list()


for(istrata in nstrata){
  
  iseed = istrata
  current_n = 10000
  current_CV = 0.15
  
  while(current_n >= 280){
    set.seed(iseed)
    
    #Create CV dataframe
    cv = list()
    for(spp in 1:ns) cv[[paste0('CV', spp)]] = current_CV
    cv[['DOM']] = 1
    cv[['domainvalue']] = 1
    cv <- as.data.frame(cv)

    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 100,
                            pops = 30,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (istrata + 1),
                            nStrata = istrata,
                            showPlot = T,
                            parallel = F)
    
    sum_stats = summaryStrata(solution$framenew,
                              solution$aggr_strata,
                              progress=FALSE) 
    
    #Update settings, res_df, and strata_list
    settings = rbind(settings,
                     data.frame(iseed = iseed,
                                nstrata = istrata,
                                cv = current_CV,
                                n = sum(sum_stats$Allocation)))
    
    strata_list = c(strata_list, sum_stats)
    res_df = cbind(res_df, solution$indices$X1)
    
    #Calculate total sample size, used to determine whether to stop 
    #optimization or repeat the optimization with a higher CV constraint
    current_n = sum(sum_stats$Allocation)
    
    #Output the results of the optimzation to the console
    print(paste0("Just Saved: ", istrata, ' strata, ',
                 current_CV*100, '% CV, ', current_n, 
                 ' sample size'))
    
    #Update the next CV level and iseed
    current_CV = current_CV + 0.005
    iseed = istrata*1000 + 1
    
    #Plot Solution
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                                 data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                              Str_no = solution$framenew$STRATO,
                                              depth = solution$framenew$X1,
                                              lon = solution$framenew$X2) )
    goa_ras = raster(goa, resolution = 5)
    goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
    plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
    
    #Save Output
    save(list = c('res_df', 'strata_list', 'settings'),
         file = paste0(github_dir, 'model_', VAST_model, '/optimization_', 
                       istrata,'_strata.RData'))
  }
  
}
