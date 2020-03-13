#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################
rm(list = ls())

library(VAST); 
library(mvtnorm); library(SamplingStrata); library(sp)
library(RColorBrewer); library(raster)

which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]

VAST_wd = c('/Users/zackoyafuso/Google Drive/VAST_Runs/',
            'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/',
            'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')[which_machine]

VAST_model = "6c"
output_wd = c(paste0('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
                     'Optimum_Allocation/model_', VAST_model),
              paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model))[which_machine]

setwd(VAST_wd)

if(!dir.exists(output_wd)) dir.create(output_wd)

load(paste0(VAST_wd, 'VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0(VAST_wd, 'VAST_output',VAST_model,'/Spatial_Settings.RData'))


load(paste0(c('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
              "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
              'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/'),
            'Extrapolation_depths.RData')[which_machine])

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
NTime = length(Years2Include)

df = df_raw = NULL

df = cbind(
  data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                          breaks = c(-171, -159, -154, -147, -140, -130), 
                          labels = c('Shumagin_1', 'Chirikof_2', 'Kodiak_3',
                                     'Yakutak_4', 'SE_5')),
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X=Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean ) )

names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth"),#, 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

for(iT in 1:NTime){
  df_raw = rbind(df_raw, cbind(
    data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                                 breaks = c(-171, -159, -154, -147, -140, -130), 
                                 labels = paste(1:5)),
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
                      X = c("depth"),#, 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.3),
                       mut_change = c(0.01, 0.1, 0.5),
                       elitism_rate = c(0.1, 0.2, 0.5),
                       dom1 = 2:5,
                       dom2 = 2:5,
                       dom3 = 2:5,
                       dom4 = 2:5,
                       dom5 = 2:5)


ns = Save$TmbData$n_c
domains = unique(df$Domain)
ndom = length(unique(frame$domainvalue))

which_runs = list(c(1:1000),
                  c(1001:5000),
                  c(5001:nrow(settings)))[which_machine]

rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid'))

# for(i in which_runs){
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
                        pops = 10,
                        elitism_rate = settings$elitism_rate[i],
                        mut_chance = settings$mut_change[i],
                        nStrata = rep(3, ndom),
                        # nStrata = unlist(settings[i, paste0('dom',1:ndom)]),
                        showPlot = F,
                        parallel = T)

strataStructure <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress=FALSE)

#   save(list=c('strataStructure', 'solution'), file = wd)
# }

