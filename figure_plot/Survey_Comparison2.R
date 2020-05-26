################################################
## 
################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(sp); library(RColorBrewer); library(raster)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]
optimization_type = c('_spatial', '_spatiotemporal')[2]
modelno = "6g"

SamplingStrata_dir = paste0(c('/Users/zackoyafuso/',
                              'C:/Users/Zack Oyafuso/',
                              'C:/Users/zack.oyafuso/')[which_machine],
                            'Downloads/SamplingStrata-master/R')

github_dir = paste0(c('/Users/zackoyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", modelno,
                   optimization_type)

###########################
## Load Current CV Simulation
###########################
load( paste0(output_wd, '/Survey_Simulation_Results.RData') )
load( paste0(output_wd, '/Stratified_RS_Simulation_Results_SurveyComp.RData') )

CV_constraints = apply(Survey_true_cv_array, MARGIN = 2:3, mean)
CV_constraints = ifelse(CV_constraints < 0.1, 0.1, CV_constraints)
stratas = c(5,10,15,20,25,30,40,50,60)
Nstrata = length(stratas)
ns = 15

##################################
## Sample Sizes
##################################
N_Matrix = matrix(nrow = Nstrata, ncol = 3)
for(isample in 1:3){
 for(istrata in 1:Nstrata){
  temp_dir = paste0(output_wd, '/Survey_Comparison/Str_', 
             stratas[istrata], 'Boat_', isample)
  load(paste0(temp_dir, '/result_list.RData'))
  N_Matrix[istrata, isample] = sum(result_list[[2]]$Allocation)
 }
}

par(mfrow = c(1,1), mar = c(5,5,1,1))
matplot(stratas, N_Matrix, type = 'l', lty = 1, col = 'black',
        ylim = c(0,1200), las = 1, xlim = c(0,65),
        xlab = 'Number of Strata', ylab = 'Total Sample Size')
matpoints(stratas, N_Matrix, pch = 16, col = 'black',
        ylim = c(0,1200), las = 1, xlim = c(0,65))
points(x = rep(65, 3), y = c(280, 550, 820), pch = 16, cex = 2)
text(x = rep(60, 3), y = c(360, 730, 850),
     paste(1:3, 'Boat'))

###############################################
## RRMSE of Estimate
###############################################

#Order of the plot
RRMSE_order = c()
for(ispp in 1:ns){
 temp_rrmse = sqrt(mean((Survey_sim_mean[1,ispp,1,] - true_mean[1,ispp])^2))/true_mean[1,ispp] 
 RRMSE_order = c(temp_rrmse, RRMSE_order)
}

{par(mfrow = c(1,4), mar = c(3,0,3,0))
 for(isample in 1:3){
  plot(1, type = 'n', xlim = c(-1,1), ylim = c(0, 15*1.5), axes = F, ann = F)
  abline(v = 0, lty = 'dashed'); box()
  mtext(side = 3, paste(isample, 'Boat'), line = 1, cex = 1.5)
  offset = 0
  axis(side = 1, at = 0, 'True Mean Density', cex.axis = 1.5)
  
  for(ispp in order(RRMSE_order, decreasing = T)){
   temp_density = density((Survey_sim_mean[1,ispp,,isample] - true_mean[1,ispp])/true_mean[1,ispp] )
   temp_density$y = offset + (temp_density$y - min(temp_density$y)) / diff(range(temp_density$y ))
   lines(temp_density$x, temp_density$y, lwd = 2)
   
   temp_density = density((STRS_sim_mean[1,ispp,3,isample,] - true_mean[1,ispp])/true_mean[1,ispp] )
   temp_density$y = offset + (temp_density$y - min(temp_density$y)) / diff(range(temp_density$y ))
   lines(temp_density$x, temp_density$y, col = 'red', lwd = 2)
   offset = offset + 1.5
  }
 }
 
 plot(1, type = 'n', xlim = c(-1,1), ylim = c(0, 15*1.5), axes = F, ann = F)
 text(x = 0, y = seq(0.25,by = 1.5, length = 15), 
      sci_names[order(RRMSE_order, decreasing = T)], font = 3, cex = 2)
 legend('top', legend = c('Survey', 'Optimization'), 
        col = c('black', 'red'), lwd = 2, lty = 1, cex = 1.5)
}

RRMSE_order = c()
for(ispp in 1:ns){
 temp_rrmse = sqrt(mean((Survey_sim_cv[1,ispp,,isample] - Survey_true_cv_array[1,ispp,isample])^2)) /Survey_true_cv_array[1,ispp,isample]
 RRMSE_order = c(temp_rrmse, RRMSE_order)
}
{par(mfrow = c(1,4), mar = c(3,0,3,0))
 for(isample in 1:3){
  plot(1, type = 'n', xlim = c(-1,1), ylim = c(0, 15*1.5), axes = F, ann = F)
  abline(v = 0, lty = 'dashed'); box()
  mtext(side = 3, paste(isample, 'Boat'), line = 1, cex = 1.5)
  axis(side = 1, at = 0, 'True CV', cex.axis = 1.5)
  offset = 0
  
  for(ispp in order(RRMSE_order, decreasing = T)){
   temp_density = density((Survey_sim_cv[1,ispp,,isample] - Survey_true_cv_array[1,ispp,isample])/Survey_true_cv_array[1,ispp,isample] )
   temp_density$y = offset + (temp_density$y - min(temp_density$y)) / diff(range(temp_density$y ))
   lines(temp_density$x, temp_density$y, lwd = 2)
   
   temp_density = density((STRS_sim_cv[1,ispp,3,isample,] - STRS_true_cv_array[1,ispp,3,isample])/STRS_true_cv_array[1,ispp,3,isample] )
   temp_density$y = offset + (temp_density$y - min(temp_density$y)) / diff(range(temp_density$y ))
   lines(temp_density$x, temp_density$y, col = 'red', lwd = 2)
   offset = offset + 1.5
  }
 }
 
 plot(1, type = 'n', xlim = c(-1,1), ylim = c(0, 15*1.5), axes = F, ann = F)
 text(x = 0, y = seq(0.25,by = 1.5, length = 15), 
      sci_names[order(RRMSE_order, decreasing = T)], font = 3, cex = 2)
 legend('top', legend = c('Survey', 'Optimization'), 
        col = c('black', 'red'), lwd = 2, lty = 1, cex = 1.5)
}

