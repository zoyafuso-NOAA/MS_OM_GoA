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
library(VAST);  library(SamplingStrata)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[1]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')
VAST_model = "6g"
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

load(paste0(github_dir, 'model_', VAST_model, '/optimization_data_model_',
            VAST_model, '.RData'))
load(paste0(github_dir, 'model_', VAST_model, '/optimization_ST_master.RData'))

frame_raw$year = rep(1:NTime, each = N)
#res_df = res_df[,-1]

strata_list = strata_list[settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
res_df = res_df[,settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
settings = subset(settings, (mut_change == 0.10 & elitism_rate == 0.10 & cv == 0.15))
settings$n = sapply(strata_list, FUN = function(x) sum(x$Allocation))

best_sol = aggregate(n ~ nstata + cv, data = settings, FUN = min)

settings = settings[row.names(settings) != '558',]
res_df = res_df[,colnames(res_df) != 'sol_558']

ids = as.numeric(rownames(res_df))
N = length(ids)
strata = c(5:20, 25, 30)
Nstrata = length(strata)
Niters = 10

#################
# true density
#################
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

sim_mean = sim_cv = array(dim = c(NTime, ns, Nstrata, Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          paste0('strata_', strata), 
                                          NULL))

for(istrata in strata) {
  print(paste('Doing Strata', istrata))
  #rownames of settings with the strata number
  row_idx = row.names(settings)[settings$nstata == istrata]
  if(length(row_idx) > 1) {
    middle_sol = floor(length(row_idx)/2)
    } else middle_sol = row_idx
  
  
  sorted_sol = order(settings[row_idx, 'n'])
  median_sol = row_idx[which(sorted_sol == middle_sol)]
  
  solno = paste0('sol_', median_sol)
  
  strata_allocation = strata_list[[solno]]$Allocation
  stratapop = strata_list[[solno]]$Population
  
  for(iyear in 1:11){
    for(iter in 1:Niters){
      sample_vec = c()
      for(i in 1:length(strata_allocation)){
        available_cells = which(res_df[,solno] == i)
        sample_cells = sample(x = available_cells, 
                              size = strata_allocation[i], 
                              replace = F)
        sample_vec = c(sample_vec, sample_cells)
      }
      sample_vec = sort(sample_vec)
      n = length(sample_vec)
      stratano =  res_df[sample_vec,solno]
      sample_df = subset(frame_raw, year == iyear)[sample_vec,]
      
      #Calculate Stratum Mean
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = mean)")
      sample_mean = eval(parse(text = stmt))[,-1]
      
      #Calculate Stratum Variance
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = var)")
      sample_var = eval(parse(text = stmt))[,-1]
      
      SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                              STATS = (stratapop/N)^2*(1 - strata_allocation/stratapop)/strata_allocation,
                              FUN = '*'))
      
      SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                               STATS = stratapop / N,
                               FUN = '*'))
      
      strata_cv = sqrt(SRS_var) / SRS_mean 
      
      sim_mean[paste0('Year_', iyear), , 
               paste0('strata_', istrata),iter] = SRS_mean
      sim_cv[paste0('Year_', iyear), , 
             paste0('strata_', istrata),iter] = strata_cv
    }
  }
}

#Save Results
save(list = c('sim_mean', 'sim_cv', 'true_mean', 'N', 'ns', 'Nstrata',
              'NTime', 'sci_names', 'strata', 'VAST_model'),
    file = paste0(github_dir, '/model_', VAST_model, '/sim_res_1.RData'))

#True CV
true_cv_array = cv_cv_array = rrmse_cv_array = 
  array(dim = c(NTime, ns, Nstrata), 
        dimnames = list(paste0('Year_', 1:NTime),
                        sci_names, 
                        paste0('strata_', strata)))

for(iyear in 1:NTime){
  for(istrata in strata){
    for(spp in sci_names){
      true_cv_array[paste0('Year_', iyear), spp, 
                    paste0('strata_', istrata)] = 
        sd(sim_mean[paste0('Year_', iyear), spp, 
                    paste0('strata_', istrata),]) / true_mean[iyear,spp]
      
      cv_cv_array[paste0('Year_', iyear), spp, 
                  paste0('strata_', istrata)] = 
        sd(sim_cv[paste0('Year_', iyear), spp, 
                  paste0('strata_', istrata), ] ) / 
        mean(sim_cv[paste0('Year_', iyear), spp, 
                    paste0('strata_', istrata), ] )
      
      rrmse_cv_array[paste0('Year_', iyear), spp, 
                     paste0('strata_', istrata)] = 
        (sum((sim_cv[paste0('Year_', iyear), spp, 
                     paste0('strata_', istrata),] - 
                true_cv_array[paste0('Year_', iyear), spp, 
                              paste0('strata_', istrata)])^2) / 
           Niters)^0.5 / mean(sim_cv[paste0('Year_', iyear), spp, 
                                     paste0('strata_', istrata) ,])
    }
  }
}

##########################
## Plot mean CV across species and years
##########################
par(mar = c(0,0,0,0), mfcol = c(5,3), oma = c(5,6,1,1))
for(spp in sci_names){ 
  plot(strata, colMeans(true_cv_array[, spp,]), type = 'n',
       ylab = '', pch = 16, axes = F, ann = F,
       las = 1, ylim = c(0,0.3))
  boxplot((true_cv_array[, spp,]), add = T, at = strata, axes = F)
  mtext(side = 3, spp, line = -2, font = 3)
  
  if(spp %in% sci_names[1:5])axis(side = 2, las = 1)
  if(spp %in% sci_names[c(5,10,15)]) axis(side = 1)
  
  box()
}
mtext(side = 2, 'True CV (boxplot distribution across years)', outer = T, line = 4)
mtext(side = 1, 'Number of Strata', outer = T, line = 3)

par(mar = c(0,0,0,0), mfcol = c(5,3), oma = c(5,6,1,1))
for(spp in sci_names){ 
  plot(strata, colMeans(apply(sim_cv[, spp,,], MARGIN = 1:2, mean)), type = 'n',
       ylab = '', pch = 16, axes = F, ann = F,
       las = 1, ylim = c(0,0.3))
  boxplot(apply(sim_cv[, spp,,], MARGIN = 1:2, mean), 
          add = T, at = strata, axes = F)
  mtext(side = 3, spp, line = -2, font = 3)
  
  if(spp %in% sci_names[1:5])axis(side = 2, las = 1)
  if(spp %in% sci_names[c(5,10,15)]) axis(side = 1)
  
  box()
}
mtext(side = 2, 'Mean CV (boxplot distribution across years)', outer = T, line = 4)
mtext(side = 1, 'Number of Strata', outer = T, line = 3)

for(spp in sci_names){ 
  plot(strata, colMeans(cv_cv_array[, spp,]), type = 'n',
       ylab = '', pch = 16, axes = F, ann = F,
       las = 1, ylim = c(0,0.6))
  boxplot((cv_cv_array[, spp,]), add = T, at = strata, axes = F)
  mtext(side = 3, spp, line = -2, font = 3)
  
  if(spp %in% sci_names[1:5]) axis(side = 2, las = 1)
  if(spp %in% sci_names[c(5,10,15)]) axis(side = 1)
  
  box()
}
mtext(side = 2, 'CV of CV (distribution across years)', outer = T, line = 4)
mtext(side = 1, 'Number of Strata', outer = T, line = 3)

for(spp in sci_names){ 
  plot(strata, colMeans(rrmse_cv_array[, spp,]), type = 'n',
       ylab = '', pch = 16, axes = F, ann = F,
       las = 1, ylim = c(0,0.6))
  boxplot((rrmse_cv_array[, spp,]), add = T, at = strata, axes = F)
  mtext(side = 3, spp, line = -2, font = 3)
  
  if(spp %in% sci_names[1:5])axis(side = 2, las = 1)
  if(spp %in% sci_names[c(5,10,15)]) axis(side = 1)
  
  box()
}
mtext(side = 2, 'RRMSE (distribution across years)', outer = T, line = 4)
mtext(side = 1, 'Number of Strata', outer = T, line = 3)

