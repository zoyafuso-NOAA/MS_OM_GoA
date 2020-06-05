#############################
## Show performance metrics across strata
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[1]
optimization_type = c('_spatial', '_spatiotemporal')[2]
VAST_model = '6g'

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", VAST_model,
                   optimization_type)

paper_dir = paste0(c('/Users/zackoyafuso/', 
                     'C:/Users/Zack Oyafuso/')[which_machine],
                   'Google Drive/MS_Optimizations/figure_plot/')
PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

load(paste0(output_wd, "/optimization_data_model_", VAST_model, ".RData"))
load(paste0(output_wd, "/Stratified_RS_Simulation_Results.RData"))
load(paste0(output_wd, "/Survey_Simulation_Results.RData"))

stratas = c(5,10,15,20,25,30,40,50,60)

##############################
## Plot
##############################
{
  png(filename = paste0(PP_dir, 'Survey_Comparison_Scheme1.png'), 
      width = 140, height = 220, units = 'mm', res = 500)
  
  layout(mat = matrix(1:(8*15), byrow = T, ncol = 8), 
         widths = c(rep(1,3),0.5,rep(1,3),1 ))
  par(mar = c(0,0,0,0), oma = c(1,3,1,0))
  
  for(ispp in 1:ns){
    
    ####################Plot True CV
    for(isample in 1:3){
      
      plot(1, type = 'n', axes = F, ann = F,xlim = c(0,12), 
           ylim = c(0,  max( STRS_true_cv_array[,ispp,,], 
                             Survey_true_cv_array[,ispp,] )) )
      if(isample == 1) {
        axis(side = 2, las = 1,
             at = pretty(x = c(0, max( STRS_true_cv_array[,ispp,,], 
                                       Survey_true_cv_array[,ispp,] )), 
                         n = 2))}
      box()
      boxplot(STRS_true_cv_array[,ispp,,isample], add = T, axes = F, 
              cex = 0.25, at = 1:9 )
      
      boxplot(Survey_true_cv_array[,ispp,isample], add = T, axes = F,cex = 0.25,
              at = 11, col = c('red', 'blue', 'grey')[isample], width = 1)
    }
    
    plot(1, type = 'n', axes = F, ann = F, 
         xlim = c(0,12), ylim = c(0, max( STRS_true_cv_array[,ispp,,])))
    
    #######################Plot RRMSE of CV
    for(isample in 1:3){
      
      plot(1, type = 'n', axes = F, ann = F,xlim = c(0,12), 
           ylim = c(0,  max( STRS_rrmse_cv_array[,ispp,,], 
                             Survey_rrmse_cv_array[,ispp,] )) )
      if(isample == 1) {
        axis(side = 2, las = 1,
             at = pretty(x = c(0,  max( STRS_rrmse_cv_array[,ispp,,], 
                                        Survey_rrmse_cv_array[,ispp,]) ),
                         n = 2))}
      box()
      boxplot(STRS_rrmse_cv_array[,ispp,,isample], add = T, axes = F,
              at = 1:9, cex = 0.25)
      
      boxplot(Survey_rrmse_cv_array[,ispp,isample], add = T, axes = F, 
              cex = 0.25,
              at = 11, col = c('red', 'blue', 'grey')[isample], width = 1)
    }

    #######################Name Plate
    plot(1, type = 'n', axes = F, ann = F, xlim = c(0,10), ylim = c(0, 10))
    text(5,5, gsub(x = sci_names[ispp], pattern = ' ', replacement = '\n' ),
         font = 3, cex = 0.8)
  }
  
  dev.off()
}
