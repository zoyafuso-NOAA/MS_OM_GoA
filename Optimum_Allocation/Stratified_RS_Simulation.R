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

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[1]
VAST_model = "6g"
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/model_', VAST_model)

VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
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
## Load data
#########################
load(paste0(github_dir, '/optimization_data_model_', VAST_model, '.RData'))
load(paste0(github_dir, '/optimization_results.RData'))

ids = as.numeric(rownames(res_df))
N = length(ids)
strata = sort(unique(settings$nstrata))
Nstrata = length(strata)
Niters = 100
sci_names = c("Atheresthes stomias", "Gadus chalcogrammus", 
              "Gadus macrocephalus", "Glyptocephalus zachirus" , 
              "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
              "Lepidopsetta bilineata", "Lepidopsetta polyxystra", 
              "Limanda aspera", "Microstomus pacificus",
              "Sebastes alutus", "Sebastes B_R", "Sebastes polyspinis", 
              "Sebastes variabilis", "Sebastolobus alascanus" )

frame_raw$year = rep(1:NTime, each = N)
res_df = res_df[,-1]

#################
# true density
#################
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

sim_mean = sim_cv = array(dim = c(NTime, ns, nrow(settings), Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          NULL, 
                                          NULL))

set.seed(134235)
for(irow in 1:nrow(settings)) {
   
   print(paste0('Settings: ', settings$nstrata[irow], ' strata, ', 
                settings$n[irow], ' samples, ', settings$cv[irow]*100, '% CV'))
   
   strata_allocation = strata_list[[irow]]$Allocation
   stratapop = strata_list[[irow]]$Population
   
   #Good strata (with > 1 cells)
   temp_strata = (1:length(strata_allocation))[strata_allocation >= 2]
   
   for(iyear in 1:NTime){
      for(iter in 1:Niters){
         
         sample_vec = c()
         for(i in temp_strata ){
            available_cells = which(res_df[,irow] == i)
            sample_cells = sample(x = available_cells, 
                                  size = strata_allocation[i], 
                                  replace = F)
            sample_vec = c(sample_vec, sample_cells)
         }
         
         sample_vec = sort(sample_vec)
         n = length(sample_vec)
         
         stratano =  res_df[sample_vec,irow]
         sample_df = subset(frame_raw, year == iyear)[sample_vec,]
         stmt = paste0('aggregate(cbind(',
                       paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                       ") ~ stratano, data = sample_df, FUN = mean)")
         sample_mean = eval(parse(text = stmt))[,-1]
         stmt = paste0('aggregate(cbind(',
                       paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                       ") ~ stratano, data = sample_df, FUN = var)")
         sample_var = eval(parse(text = stmt))[,-1]
         
         temp_strata_allocation = strata_allocation[temp_strata]
         temp_stratapop = stratapop[temp_strata]
         
         SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                                 STATS = (temp_stratapop/N)^2*(1 - temp_strata_allocation/temp_stratapop)/temp_strata_allocation,
                                 FUN = '*'))
         
         SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                                  STATS = temp_stratapop / N,
                                  FUN = '*'))
         
         strata_cv = sqrt(SRS_var) / SRS_mean 
         
         sim_mean[paste0('Year_', iyear), , irow,iter] = SRS_mean
         sim_cv[paste0('Year_', iyear), , irow,iter] = strata_cv
         
      }
   }
}

#################################
## Simulation Metrics
#################################
#True CV
true_cv_array = cv_cv_array = rrmse_cv_array = 
   array(dim = c(NTime, ns, nrow(settings)), 
         dimnames = list(paste0('Year_', 1:NTime), sci_names, NULL ))

for(iyear in 1:NTime){
   for(irow in 1:nrow(settings) ){
      for(spp in sci_names){
         true_cv_array[paste0('Year_', iyear), spp, 
                       irow] = 
            sd(sim_mean[paste0('Year_', iyear), spp, 
                        irow,]) / true_mean[iyear,spp]
         
         cv_cv_array[paste0('Year_', iyear), spp, 
                     irow] = 
            sd(sim_cv[paste0('Year_', iyear), spp, irow, ] ) / 
            mean(sim_cv[paste0('Year_', iyear), spp, irow, ] )
         
         rrmse_cv_array[paste0('Year_', iyear), spp, irow] = 
            (sum((sim_cv[paste0('Year_', iyear), spp, irow,] - 
                     true_cv_array[paste0('Year_', iyear), spp, irow])^2) / Niters)^0.5 / mean(sim_cv[paste0('Year_', iyear), spp, irow ,])
      }
   }
}

#######################
## Save results
#######################
for(ivar in  c('cv_cv_array', 'rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv')){
   assign(x=paste0('STRS_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, '/Stratified_RS_Simulation_Results.RData'),
     list = c(paste0('STRS_', c('cv_cv_array', 'rrmse_cv_array', 
                               'true_cv_array', 'sim_mean', 'sim_cv')),
              'true_mean', 'sci_names', 'NTime', 'ns', 'Niters', 'N'))


