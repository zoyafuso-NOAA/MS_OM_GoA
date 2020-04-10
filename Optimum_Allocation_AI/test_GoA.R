##############################
## Bethel Algorithm on GoA current strata
##############################
rm(list = ls())

############################
## Import Libraries
############################
library(rgdal); library(raster); library(rgeos); library(tidyr)
library(SamplingStrata)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI' = 3)[3]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine], 'GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/VAST_Runs/',
                    'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/',
                    'C:\\Users\\zack.oyafuso\\Desktop\\VAST_Runs\\'),
                  'VAST_output', modelno)[which_machine]

############################
## Load Data
############################
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
load(paste0(github_dir, '/Optimum_Allocation/model_', modelno, 
            '/optimization_data_model_', modelno, '.RData'))

#Survey data
survey_data = read.csv(paste0(github_dir, '/data/data/',
                              'GOA_multspp_with_strata.csv'))
survey_data = subset(survey_data, SPECIES_NAME != 'Anoplopoma fimbria')
survey_data$SPECIES_NAME = droplevels(survey_data$SPECIES_NAME)

## Constants
sci_names = levels(Data_Geostat$spp)
N = nrow(frame)

#Calculate weights of each stratum
stratapop = table(gulf_of_alaska_grid$GOA_STRATUM)

#Calculate samples allocated across strata across years
samples_by_str = with(subset(survey_data, 
                             SPECIES_NAME == 'Sebastolobus alascanus'),
                      table(YEAR, STRATUM))
strata = colnames(samples_by_str)


######

CVs = seq(0.10, 0.35, by = 0.01)
Years = paste(sort(unique(survey_data$YEAR)))
sample_allocation = array(data = 0, dim = c(NTime, length(CVs), length(strata)),
                          dimnames = list(Years, CVs, strata))

for(iyear in Years){
  sample_df = subset(survey_data, YEAR == iyear)
  temp_strata = paste(sort(unique(sample_df$STRATUM)))
  sample_mean = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM + YEAR,
                                        data = sample_df, FUN = mean, drop = F),
                       key = SPECIES_NAME, value = CPUE)[,-(1:2)]
  names(sample_mean) = paste0('M', 1:ns)
  
  sample_var = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM + YEAR,
                                       data = sample_df, FUN = var, drop = F),
                      key = SPECIES_NAME, value = CPUE)[,-(1:2)]
  names(sample_var) = paste0('S', 1:ns)
  
  sample_mean[is.na(sample_var)] = 0
  sample_var[is.na(sample_var)] = 0
  temp_stratapop = stratapop[temp_strata]
  
  
  df = cbind(data.frame(stratum = temp_strata,
                        N = as.vector(temp_stratapop)),#,
             #X1 = factor(1:length(temp_strata))),
             sample_mean, sqrt(sample_var),
             data.frame(cens = 0,
                        cost = 1,
                        DOM1 = 'tot'))
  
  for(icv in CVs){
    stmt = paste0('cbind(', paste0("CV", 1:ns, '=', icv, collapse = ', '), ')' )
    CV = eval(parse(text = stmt))
    errors = cbind(data.frame(DOM = 'DOM1'), CV, domainvalue = 1)
    
    n = bethel(stratif = df, errors = errors, printa=TRUE, epsilon = 1e-11, maxiter = 200)
    sample_allocation[iyear, paste(icv), temp_strata] = as.numeric(n)
  }
  
}




total_sample_size = apply(sample_allocation, MARGIN = 2, sum)
plot(x = CVs, y = total_sample_size, pch = 16, las = 1, ylim = c(0, 2000))
lines(x = CVs, y = total_sample_size)

attributes(n)$outcv

