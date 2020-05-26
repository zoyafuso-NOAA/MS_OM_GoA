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
                   optimization_type, '/Flexible_Optimization/')

stratas = c(5,10,15,20,25,30,40,50,60)
NStrata = length(stratas)
ns = 15
spp_cv = samplesizes = list()

for(istrata in 1:2){
 temp_strata = paste0('Str_', stratas[istrata])
 runs = grep(x = dir(output_wd, full.names = T), 
             pattern = paste0('Str_', stratas[istrata]),
             value = T)
 
 nruns = length(runs)
 for(irun in 1:nruns){
  load( paste0(runs[irun], '/result_list.RData') )
  samplesizes[[temp_strata]]$n = c(samplesizes[[temp_strata]]$n, result_list$n) 
  spp_cv[[temp_strata]]$cv = rbind(spp_cv[[temp_strata]]$cv, result_list[[3]])
 }
}

par(mar = c(5,5,1,1))
matplot( t(spp_cv[['Str_10']]$cv[,order(spp_cv[['Str_10']]$cv[1,])] ), 
         type = 'b', lty = 1, pch = paste(1:nruns),
         las = 1, xlab = 'Species', ylim = c(0,0.3),
         ylab = 'Expected Spatiotemporal CV',
         lwd = c(1,3,1,1,3,1), cex = c(1,1.5,1,1,1.5,1),
         col = c('darkgrey','black', 'darkgrey', 'darkgrey', 'black',
                 'darkgrey') )
         
plot(samplesizes[['Str_10']]$n, pch = paste(1:nruns), type = 'b', cex = 2,
     xlab = 'Run Number', ylab = 'Total Sample Size', las = 1, ylim = c(0,850))
abline(h = c(280, 550, 820), col = 'darkgrey', lty = 'dashed')
