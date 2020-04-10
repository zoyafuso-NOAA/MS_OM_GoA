#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine],
                    'GitHub/MS_OM_GoA/Optimum_Allocation/', 'model_', modelno, '/')

load(paste0(github_dir, 'Survey_Simulation_Results.RData'))
load(paste0(github_dir, 'Simple_RS_Simulation_Results.RData'))
load(paste0(github_dir, 'Stratified_RS_Simulation_Results.RData'))
load(paste0(github_dir, 'optimization_results.RData'))


plot_set = data.frame(spp = sci_names,
                      true_ymax = c(0.10, 0.25, 0.10, 0.10, 0.15,
                                    0.06, 0.15, 0.25, 0.50, 0.20,
                                    0.25, 0.20, 0.45, 0.40, 0.20),
                      rrmse_ymax =c(0.20, 0.45, 0.35, 0.35, 0.25,
                                    0.25, 0.40, 0.35, 0.30, 0.50,
                                    0.33, 0.40, 0.75, 0.80, 0.50))

par(mfrow = c(5,2), mar = c(0,5,3,0), oma = c(2,0,0,0))
for(spp in 1:15){
  for(imetric in c('true', 'rrmse')){
    survey_result = get(paste0('survey_', imetric, '_cv_array'))
    SRS_result = get(paste0('SRS_', imetric, '_cv_array'))
    STRS_result = get(paste0('STRS_', imetric, '_cv_array'))
    
    #Empty plot for True CV
    plot(1, type = 'n', xlim = c(0,17), 
         ylim = c(0, plot_set[spp, paste0(imetric, '_ymax')] ), 
         ann = F, axes = F)
    box()
    mtext(side = 3,  paste(spp, sci_names[spp]))
    #axis(side = 1, at = c(2:3, 5, 7), labels = NA)
    axis(side = 2, las = 1)
    
    #Simulated Survey Strata
    boxplot( survey_result[,spp], add = T, width = 1, axes = F, at = 1)
    
    #Simulated Simple Random Sampling, 350 and 550 samples
    boxplot( SRS_result[,spp,c('size_550', 'size_350')], add = T, 
             axes = F, at = 3:4)
    
    #Simulated Stratified Random Sampling, optimized strata
    which_strata = which(settings$nstrata == 5)
    row_idx = which_strata[which.min(settings[which_strata,'n'])]
    
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 6)
    
    which_strata = which(settings$nstrata == 20)
    row_idx = which_strata[which.min(settings[which_strata,'n'])]
    
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 12)
    
    which_strata = which(settings$nstrata == 40)
    row_idx = which_strata[which.min(settings[which_strata,'n'])]
    
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 15)
    
  }
}

