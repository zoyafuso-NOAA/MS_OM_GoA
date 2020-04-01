######################################
## Calculate sampling CV
######################################

#####################################
## Optimal Solutions from 5-20 strata
#####################################
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
frame_raw$year = rep(1:11, each = nrow(frame))

ns = Save$TmbData$n_c
sci_names = Save$Spp


rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid', 'iT', 'df_raw'))


load('optimization_ST_master.RData')

strata_list = strata_list[settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
res_df = res_df[,settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
settings = subset(settings, (mut_change == 0.10 & elitism_rate == 0.10 & cv == 0.15))
settings$n = sapply(strata_list, FUN = function(x) sum(x$Allocation))

best_sol = aggregate(n ~ nstata + cv, data = settings, FUN = min)




ids = as.numeric(rownames(res_df))
N = length(ids)

cv_array = array(dim = c(11, ns, 16, 100), 
                 dimnames = list(paste0('Year_', 1:11),
                                 sci_names, paste0('strata_', 5:20), NULL))


for(istrata in 1:nrow(best_sol)) {
  #rownames of settings with the strata number
  row_idx = row.names(settings)[settings$nstata == best_sol$nstata[istrata]]
  best_sol_idx = which.min(settings[row_idx, 'n'])
  idx = row_idx[best_sol_idx]
  solno = paste0('sol_', idx)
  
  strata_allocation = strata_list[[solno]]$Allocation
  stratapop = strata_list[[solno]]$Population
  for(iyear in 1:11){
    for(iter in 1:100){
      sample_vec = c()
      for(i in 1:length(strata_allocation)){
        available_cells = which(res_df[,solno] == i)
        sample_cells = sample(x = available_cells, 
                              size = strata_allocation[i], 
                              replace = F)
        sample_vec = c(sample_vec, sample_cells)
      }
      
      sample_vec = sort(sample_vec)
      stratano =  res_df[sample_vec,solno]
      sample_df = subset(frame_raw, year == iyear)[sample_vec,]
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = mean)")
      sample_mean = eval(parse(text = stmt))[,-1]
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = var)")
      sample_var = eval(parse(text = stmt))[,-1]
      
      SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                              STATS = 1/strata_allocation * (stratapop - strata_allocation)/stratapop * (stratapop/N)^2,
                              FUN = '*'))

      SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                               STATS = stratapop / N,
                               FUN = '*'))
      
      strata_cv = sqrt(SRS_var) / SRS_mean * 100
      
      cv_array[paste0('Year_', iyear), , 
               paste0('strata_', (5:20)[istrata]),iter] = strata_cv
      
    }
  }
}

for(ispp in sci_names){
  plot(1, type = 'n', xlim = c(5,20), ylim = c(0,5), las = 1, main = ispp)
  
  for(istrata in 5:20){
    boxplot(at = istrata, apply(cv_array[,ispp,paste0('strata_', istrata),], MARGIN = 2, FUN = function(x)sd(x) ), add = T, axes = F )
  }
}

for(ispp in sci_names){
  plot(1, type = 'n', xlim = c(5,20), ylim = c(0,20), las = 1, main = ispp)
  
  for(istrata in 5:20){
    boxplot(at = istrata, apply(cv_array[,ispp,paste0('strata_', istrata),], MARGIN = 2, FUN = mean), add = T, axes = F )
  }
}

 
100*apply(cv_array[,,paste0('strata_', istrata),], MARGIN = 1:2, FUN = sd) / 
  apply(cv_array[,,paste0('strata_', istrata),], MARGIN = 1:2, FUN = mean)

boxplot(t(cv_array[,'Atheresthes stomias',paste0('strata_', istrata),]))
