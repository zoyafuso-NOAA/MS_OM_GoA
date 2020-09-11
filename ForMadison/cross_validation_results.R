###############################################################################
## Project:       Cross-Validation Results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For each CV run, calculate relative root mean squre error of
##                density predictions
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(VAST)
library(RANN)
library(tidyr)

##################################################
####   Set up directores
##################################################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[3]

VAST_dir = c("",
             'C:/Users/Zack Oyafuso/Google Drive/GOA_VAST_Runs/Single_Species/',
             "C:/Users/zack.oyafuso/Desktop/VAST_Runs/Single_Species/")[which_machine]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/', 
                      'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                    "GitHub/Optimal_Allocation_GoA/")

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'model_11/optimization_data.RData'))

which_spp = c( 
  'Atheresthes stomias', 
  'Gadus chalcogrammus', 
  'Gadus macrocephalus',
  
  'Glyptocephalus zachirus', 
  'Hippoglossoides elassodon', 
  'Hippoglossus stenolepis', 
  
  'Lepidopsetta bilineata', 
  'Lepidopsetta polyxystra',
  'Microstomus pacificus',
  
  'Sebastes alutus',
  'Sebastes B_R',
  'Sebastes brevispinis',
  
  'Sebastes polyspinis', 
  'Sebastes variabilis',
  'Sebastolobus alascanus')[]
ns = length(which_spp)

GOA_grid = make_extrapolation_info(Region = 'Gulf_of_Alaska')$Data_Extrap
N = nrow(GOA_grid)

##################################################
####   Result Objects
##################################################
NFold = 10
NTime = 11
CV_df = expand.grid(species = which_spp, 
                    depth = c(T,F),
                    fold = 1:NFold,
                    stringsAsFactors = F)

CV_df[,c('max_grad', 'pdHess', 'bound_check', 'pred_nll', paste0("year", 1:NTime))] = NA

for(irow in 1:nrow(CV_df)){
  
  #Load fitted object
  result_dir = paste0(VAST_dir, CV_df$species[irow], 
                      ifelse(CV_df$depth[irow],  '_depth', ''), '/')
  filename = paste0(result_dir, "CV_", CV_df$fold[irow], '/fit.RData')
  
  if ( file.exists(filename) ){
    load(filename)
    
    #Final Gradient
    CV_df$max_grad[irow] = 
      max(abs(fit_new$parameter_estimates$diagnostics$final_gradient))
    
    #Check whether hessian matrix is positive definite
    CV_df$pdHess[irow] = fit_new$parameter_estimates$SD$pdHess
    
    #check_fit chekcs bounds, TRUE is bad and FALSE is good
    CV_df$bound_check[irow] = (check_fit(fit_new$parameter_estimates))
    
    CV_df$pred_nll[irow] = fit_new$Report$pred_jnll
    
    #Extract incides of withheld data
    withheld_idx = which(fit_new$data_list$PredTF_i == T)
    
    #Withheld data locations, species/year indices, observed CPUE
    withheld_df = data.frame(
      idx =  withheld_idx,
      E_km = fit_new$spatial_list$loc_i[withheld_idx,'E_km'],
      N_km = fit_new$spatial_list$loc_i[withheld_idx,'N_km'],
      year = 1+fit_new$data_list$t_iz[withheld_idx,1],
      obs_density = (fit_new$data_frame$b_i/fit_new$data_frame$a_i)[withheld_idx]
    )
    
    #Extrapolation Grid locations
    loc_g = fit_new$spatial_list$loc_g
    
    #which extrapoaltion grid cells are closest to each withheld data location 
    grid_idx = RANN::nn2(query = withheld_df[,c('E_km', 'N_km')], 
                         data = loc_g, 
                         k = 1)$nn.idx
    
    for(jrow in 1:nrow(withheld_df)) {
      temp_density = fit_new$Report$D_gcy[grid_idx[jrow], ,withheld_df$year[jrow]]
      withheld_df$pred_density[jrow] =  temp_density
    }
    
    #Calculate mean predicted density for calculation of RRMSE
    mean_pred_density = aggregate(obs_density ~ year, data = withheld_df, 
                                  FUN = mean)
    
    for(itime in 1:NTime){
      split_df = subset(withheld_df, year ==  unique(withheld_df$year)[itime])
      temp_RMSE = sqrt(mean((split_df$obs_density - split_df$pred_density)^2))
      temp_mean_pred_density = mean_pred_density[itime, 'obs_density']
      CV_df[irow,paste0('year',itime)] = temp_RMSE / temp_mean_pred_density
    }
    
    print(paste0('Done with: ', CV_df$species[irow], ', ', 
                 ifelse(CV_df$depth[irow], 'Depth, ', 'No Depth, '),
                 'Fold Number ', CV_df$fold[irow]))
  }
}

##################################################
####   Calculate Mean relatie root mean square error of predictions
####   for models that include and don't include depth
##################################################
CV_df$RRMSE = rowMeans(CV_df[,paste0('year', 1:NTime)])

RRMSE = tidyr::spread(data = aggregate(RRMSE ~ species + depth,
                                       data = CV_df, 
                                       FUN = mean), 
                      key = 'depth',  
                      value = 'RRMSE')
names(RRMSE)[-1] = c('No_Depth', 'Depth')

##################################################
#### 
##################################################
tidyr::spread(data = aggregate(pred_nll ~ species + depth,
                               data = CV_df,
                               FUN = sum),
              key = "depth",
              value = "pred_nll")

##################################################
####   Create the result object that would go into the optimizations
##################################################
D_gcy = array(dim = c(N, ns, 24), dimnames = list(NULL, which_spp, NULL))

for(ispp in 1:ns){
  which_model_idx = which.min(RRMSE[ispp, -1])
  
  result_dir = paste0(VAST_dir, RRMSE$species[ispp], 
                      c('', '_depth')[which_model_idx], '/')
  load(paste0(result_dir, '/fit.RData'))
  
  D_gcy[,ispp,] = fit$Report$D_gcy[,1,]
}

rm(fit)
fit = list(Report = list(D_gcy = D_gcy))

##################################################
####   Save
##################################################
if(!dir.exists(paste0(dirname(VAST_dir), '/VAST_output11'))) 
  dir.create(paste0(dirname(VAST_dir), '/VAST_output11'))

save('fit', file = paste0(dirname(VAST_dir), '/VAST_output11/fit.RData') )
