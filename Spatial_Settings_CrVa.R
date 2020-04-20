###############################
## Spatial Settings for VAST
###############################
rm(list = ls())

library(VAST)

###################################
## Set up directories
###################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2, 'VM' = 3)[2]

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output7/')

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

which_spp = unlist(spp_df['6j',])

Data_Geostat = subset(Data_Geostat, spp %in% names(which_spp)[which_spp])
Data_Geostat$spp = droplevels(Data_Geostat$spp)

################################
## Assign 10 fold partitions of the data
################################
# Generate partitions in data
n_fold = 10

#Sort Data_Geostat
Data_Geostat = Data_Geostat[order(Data_Geostat$Year, Data_Geostat$spp),]
Data_Geostat$id = 1:nrow(Data_Geostat)

set.seed(2342)
foldno = lapply(X = split.data.frame(Data_Geostat, f = Data_Geostat$Year),
                FUN = function(test) {
                  row.idx = matrix(data = test$id, ncol = 15)
                  fold_no = sample(x = 1:n_fold, 
                                   size = nrow(row.idx), 
                                   replace = T)
                  return(split(row.idx, fold_no))
                })

test = lapply(X = foldno,
              FUN = function(test){
                lapply(test, FUN = function(x) Data_Geostat$id[x])
              })

for(iyear in names(test)){
  for(ifold in paste(1:10)){
    Data_Geostat[test[[iyear]][[ifold]],'fold'] = as.integer(ifold)
  }
}

save(list = c('Data_Geostat', 'n_fold'),
     file = paste0(VAST_dir,'/Spatial_Settings_CrVa.RData') )
