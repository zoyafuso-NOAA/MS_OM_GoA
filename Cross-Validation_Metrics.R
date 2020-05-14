################################
## Cross Validation Metrics
################################

library(VAST); library(RANN)

################################
## Set up directories
################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2)[2]
modelno = '8a'
github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', modelno, '/')

################################
## Import Data
## Grid locations
################################
load(paste0(dirname(VAST_dir), '/Spatial_Settings_CrVa.RData'))


################################
##
################################
n_fold = 3
RRMSE = vector(length = n_fold)
SqErs = list()

for(ifold in seq(n_fold)){
  load( paste0(VAST_dir, 'CV_', ifold, '/fit.RData') )
  
  #Extract the indices for year and species
  new_years = sapply(Data_Geostat[Data_Geostat$fold == ifold, "Year"], 
                     function(x) which(x == fit_new$year_labels))
  new_spp = sapply(Data_Geostat[Data_Geostat$fold == ifold, "spp"], 
                   function(x) which(x == levels(Data_Geostat$spp)))
  
  #Calculate the index nearest grid cell for each datum
  new_locs = Data_Geostat[Data_Geostat$fold == ifold, c('Lon', 'Lat')]
  new_locs_UTM = project_coordinates(X = new_locs[, "Lon"], 
                                     Y = new_locs[, "Lat"], 
                                     flip_around_dateline = FALSE)
  grid_idx = as.vector(RANN::nn2(query=new_locs_UTM,
                                 data = fit_new$spatial_list$loc_g,
                                 k = 1)$nn.idx)
  
  #Extract predicted density
  new_density = vector(length = length(grid_idx))
  for(i in 1:length(new_density)){
    new_density[i] = fit_new$Report$D_gcy[grid_idx[i], new_spp[i], new_years[i]]
  }
  
  #Extract observed density
  new_obs_density = with(subset(Data_Geostat,fold == ifold),
                         Catch_KG/AreaSwept_km2 )
  
  #Calculate relative root mean square error
  RRMSE[ifold] = sqrt(mean((new_density-new_obs_density)^2)) / mean(new_density)
  SqErs[[ifold]] = (new_density-new_obs_density)^2
}

mean(RRMSE)

ifold = 1
hist(SqErs[[ifold]])



