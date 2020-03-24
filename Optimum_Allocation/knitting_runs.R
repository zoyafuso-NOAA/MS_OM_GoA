###########################
## Thread spatiotemporal results back together
###########################
rm(list = ls())
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
modelno = '6g'
results_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/')[which_machine],
                     'GitHub/MS_OM_GoA/Optimum_Allocation/', 'model_',
                     modelno)
setwd(results_dir)

res_files = dir(pattern = 'optimization_spatiotemporal_')

load('optimization.RData')
master_res = as.data.frame(res_df)[,1:2]
master_strata_list = list()

for(ifile in res_files){
 load(ifile)
 
 filenos = strsplit(x = gsub(gsub(ifile, pattern = 'optimization_spatiotemporal_', replacement = ''), pattern = '.RData', replacement = ''), split = '-')
 
 filenos = sort(unlist(lapply(filenos, FUN = function(x) x[1]:x[2])))
 
 res_df = as.data.frame(res_df)[,-1:-2]
 names(res_df) = paste0('sol_', filenos)
 
 master_res = cbind(master_res, res_df)
 master_strata_list[filenos] = strata_list[filenos]
}

master_res = master_res[,-c(1:2)]
names(master_strata_list) = paste0('sol_', 1:length(master_strata_list))
settings = settings[(1:length(master_strata_list))[!sapply(master_strata_list, is.null)],]
master_strata_list = master_strata_list[!sapply(master_strata_list, is.null)]

res_df = master_res
strata_list = master_strata_list

save(list = c('strata_list', 'res_df', 'settings', 
              'frame', 'ns', 'NTime', 'VAST_model'), 
     file = paste0(results_dir, '/optimization_ST_master.RData') )

