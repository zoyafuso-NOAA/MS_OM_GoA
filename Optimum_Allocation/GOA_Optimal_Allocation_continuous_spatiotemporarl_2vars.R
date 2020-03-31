#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################
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
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]

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
VAST_dir = paste0(c('', '', 
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

#########################
## Load VAST products
#########################
load(paste0(VAST_dir, '/VAST_MS_GoA_Run.RData'))
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
load(paste0(dirname(github_dir), '/Extrapolation_depths.RData'))

#########################
## Index years that had data
#########################
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
NTime = length(Years2Include)

##########################
## Create the data inputs to SamplingStrata
##########################
df = df_raw = NULL

df = cbind(
  data.frame(Domain = 1,
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X=Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean ) )

names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth", 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

for(iT in 1:NTime){
  df_raw = rbind(df_raw, cbind(
    data.frame(Domain = 1,
               x = 1:Save$TmbData$n_g,
               year = iT,
               lat = Extrapolation_depths$N_km,
               lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
               depth = Extrapolation_depths$depth),
    Save$Report$Index_gcyl[,,Years2Include[iT],] )
  )
}
names(df_raw)[-(1:6)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')
frame_raw <- buildFrameDF(df = df_raw,
                          id = "x",
                          X = c("depth", 'lon'),
                          Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                          domainvalue = "Domain")

############################
## Settings for optimizer
############################
settings = rbind(expand.grid(cv = c(0.2, 0.15),
                             mut_change = c(0.1, 0.01),
                             elitism_rate = c(0.2, 0.1),
                             nstata = c(5,7,10),
                             iter = 1:10),
                 expand.grid(cv = c(0.2, 0.15),
                             mut_change = 0.1,
                             elitism_rate = 0.1,
                             nstata = c(6,8,11:20),
                             iter = 1:10),
                 expand.grid(cv = c(0.2, 0.15),
                             mut_change = 0.1,
                             elitism_rate = 0.1,
                             nstata = 9,
                             iter = 1:10),
                 expand.grid(cv = c(0.2, 0.15),
                             mut_change = 0.1,
                             elitism_rate = 0.1,
                             nstata = 5:20,
                             iter = 11:20)
                 
)

ns = Save$TmbData$n_c

rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid', 'ifile', 'iT', 'Extrapolation_depths', 'df_raw'))

res_df = as.matrix(frame[,c('id', 'domainvalue')])
strata_list = list()

iter_range = unlist(list('Zack_MAC'= NA, 'Zack_PC' = 241:320,
                         'Zack_GI_PC'=321:400, 'VM' = 401:480)[which_machine])

for(ii in iter_range){
  
  cv = list()
  for(spp in 1:ns) cv[[paste0('CV', spp)]] = settings$cv[ii]
  cv[['DOM']] = 1
  cv[['domainvalue']] = 1
  cv <- as.data.frame(cv)
  
  set.seed(1234 + ii)
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 50,
                          pops = 30,
                          elitism_rate = settings$elitism_rate[ii],
                          mut_chance = settings$mut_change[ii],
                          nStrata = settings$nstata[ii],
                          showPlot = T,
                          parallel = F)
  
  strata_list[[ii]] =  summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress=FALSE) 
  
  res_df = cbind(res_df, solution$framenew$STRATO)
  
  save(list = c('strata_list', 'res_df', 'settings', 
                'frame', 'ns', 'NTime', 'VAST_model'), 
       file = paste0(output_wd, '/optimization_spatiotemporal_',
                     min(iter_range), '-', ii,'.RData') )
  
  #Plot
  goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                               data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                            Str_no = solution$framenew$STRATO,
                                            depth = solution$framenew$X1,
                                            lon = solution$framenew$X2) )
  goa_ras = raster(goa, resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
  plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
}

