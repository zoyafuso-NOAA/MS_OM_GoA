library(reshape2)

EBS = read.csv('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/data-raw/cpue_EBSshelf_all_spp.csv')

spp = c(10110, 10130, 10210)

df = subset(EBS,
             YEAR > 1986 
             &!(STATIONID %in% 
                  c('AZ0504', 'GF1918', 'GF2019', 'GF2120', 'GF2221', 'HG1918', 
                    'HG2019', 'HG2120','HG2221', 'IH1918', 'IH2019', 'IH2120', 
                    'IH2221', 'ON2524', 'ON2625', 'PO2423', 'PO2524', 'PO2625',
                    'PO2726', 'QP2423', 'QP2524', 'QP2625', 'QP2726', 'JI1918',
                    'JI2019', 'JI2120', 'JI2221'))
            & SPECIES_CODE %in% spp,  
             select = c('YEAR', 'STATIONID', 'SPECIES_CODE', 'WGTCPUE'))

df2 = data.frame()

 for(id in sort(unique(df$STATIONID)[unique(df$STATIONID) != 'J-13'] )){
#for(id in unique(df$STATIONID) ){
  temp_df = subset(df, STATIONID == id)
  temp_df = temp_df[!duplicated(temp_df),]
  
  temp_df = spread(temp_df, key = SPECIES_CODE, value = WGTCPUE)
  temp_df = apply(X = temp_df[,paste(spp)],
                  MARGIN = 2, 
                  FUN = function(x) (x - min(x))/diff(range(x)) )
  temp_df[is.nan(temp_df)] = 0
 
  spp_covar = cov(temp_df[,paste(spp)])
  
  df2 = rbind(df2,
              data.frame(
                STATIONID = id,
                tot_var =   sum( spp_covar[upper.tri(x = spp_covar, diag = TRUE)] ),
                tot_mean = sum(x = colMeans(temp_df[,paste(spp)])) )
              ) 


}

stations = read.csv('10110_stations_.csv')
stations = stations[stations$STATIONID != 'J-13',]

df2[,c('lat', 'long')] = stations[,c('lat', 'long')]
plot(lat ~ long, data = df2, cex = tot_var^2 * 50)
plot(lat ~ long, data = df2, cex = tot_mean^2 * 10)

write.csv(x = df2, file = 'optim_data.csv', row.names = F)
