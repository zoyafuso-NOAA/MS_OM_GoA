# Set local working directory (change for your machine)
load("C:/Users/zack.oyafuso/Desktop/VAST_Runs/VAST_output6g/VAST_MS_GoA_Run.RData")
load("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/VAST_output6j/Model_Settings.RData")
load("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/VAST_output6j/Spatial_Settings.RData")

# Load packages
library(TMB)
library(VAST)

# Generate partitions in data
n_fold = 10

#Sort Data_Geostat
Data_Geostat$id = 1:nrow(Data_Geostat)


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

lapply(split(Data_Geostat, Data_Geostat$Year), FUN = function(x) table(x$fold, x$spp))


###
# Generate partitions in data
n_fold = 10
Partition_i = sample( 1:n_fold, size=nrow(example$sampling_data), replace=TRUE )
prednll_f = rep(NA, n_fold )

# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Partition_i==fI, TRUE, FALSE )
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model( "settings"=settings, "Lat_i"=example$sampling_data[,'Lat'],
                       "Lon_i"=example$sampling_data[,'Lon'], "t_i"=example$sampling_data[,'Year'],
                       "c_i"=rep(0,nrow(example$sampling_data)), "b_i"=example$sampling_data[,'Catch_KG'],
                       "a_i"=example$sampling_data[,'AreaSwept_km2'], "v_i"=example$sampling_data[,'Vessel'],
                       "PredTF_i"=PredTF_i, "Parameters"=ParHat, "getsd"=FALSE )
  
  # Save fit to out-of-bag data
  prednll_f[fI] = fit_new$Report$pred_jnll
}

# Check fit to all out=of-bag data and use as metric of out-of-bag performance
sum( prednll_f )