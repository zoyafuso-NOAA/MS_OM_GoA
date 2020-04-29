###############################
## Spatial Settings for VAST
###############################
rm(list = ls())

library(VAST)

###################################
## Set up directories
###################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2, 'VM' = 3)[1]

modelno = '8'

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', modelno, '/')

## Import Data
data = read.csv(file = paste0(github_dir, 'data/data/GOA_multspp.csv') )

# Prepare the Data-frame for catch-rate data
Data_Geostat = data.frame( "spp"=data$SPECIES_NAME,
                           "Year"=data$YEAR,
                           "Catch_KG"=data$WEIGHT,
                           "AreaSwept_km2"=data$EFFORT,
                           "Vessel"=0,
                           "Lat"=data$LATITUDE,
                           "Lon"=data$LONGITUDE,
                           "LOG_DEPTH" = data$DEPTH, #centered log_depth
                           "LOG_DEPTH2" = data$DEPTH^2 )
rm(data)

#Drop factor levels of unused Species
spp_df = read.csv(paste0(github_dir, "spp_df.csv"), check.names=F, header = T, 
                  row.names = 'modelno')

which_spp = unlist(spp_df[modelno,])

Data_Geostat = subset(Data_Geostat, spp %in% names(which_spp)[which_spp])
Data_Geostat$spp = droplevels(Data_Geostat$spp)

################################
## Assign 10 fold partitions of the data
################################
# Generate partitions in data
n_fold = 10
years = paste0(unique(Data_Geostat$Year))
NTime = length(unique(Data_Geostat$Year))
ns = length(unique(Data_Geostat$spp))

#Sort Data_Geostat
Data_Geostat = Data_Geostat[order(Data_Geostat$Year, Data_Geostat$spp),]
Data_Geostat$latlon = paste0(Data_Geostat$Lat, Data_Geostat$Lon)

set.seed(2342)
foldno = lapply(X = split.data.frame(Data_Geostat, f = Data_Geostat$Year),
                FUN = function(test) {
                  unique_loc = unique(test$latlon)
                  fold_no = sample(x = 1:n_fold, 
                                   size = length(unique_loc), 
                                   replace = T)
                  return(split(unique_loc, fold_no))
                })



for(iyear in years){
  for(ifold in paste(1:n_fold)){
    Data_Geostat[Data_Geostat$latlon %in% foldno[[iyear]][[ifold]] ,
                 'fold'] = as.integer(ifold) 
  }
}

save(list = c('Data_Geostat'),
     file = paste0(VAST_dir,'/Spatial_Settings_CrVa.RData') )
