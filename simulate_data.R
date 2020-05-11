###############################
## Simulate Data
###############################
if(!dir.exists(paste0(VAST_dir, 'VAST_output', modelno, '/simulated_data'))){
  
  for(datatype in c('simulated_data', 'reduced_simulated_data')){
    dir.create(paste0(VAST_dir, 'VAST_output', modelno,'/', datatype))
    for(itype in 1:3){
      dir.create(paste0(VAST_dir, 'VAST_output', modelno, '/', datatype,
                        '/type', itype, '/'))
    }
  }
}
Niters = 10
for(i in 1:Niters){
  for(itype in 1:3){
    temp_sim = simulate_data(fit, type = itype) 
    
    save('temp_sim', 
         file = paste0(VAST_dir, 'VAST_output', modelno, 
                       '/simulated_data/type', itype, '/sim_', i, '.RData'))
    
    temp_sim_survey = temp_sim[[1]][c('Index_gcyl', 'D_gcy')]
    save('temp_sim_survey', 
         file = paste0(VAST_dir, 'VAST_output', modelno, 
                       '/reduced_simulated_data/type',itype,'/sim_',i,'.RData'))
  }
}

rm(Niters, i, itype, temp_sim, temp_sim_survey)
