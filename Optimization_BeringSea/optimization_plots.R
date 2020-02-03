#################################
## Optimization Plots
################################

setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/Optimization_BeringSea')

load('optimization_results.RData')
stations = read.csv('10110_stations_.csv')
stations = stations[stations$STATIONID != 'J-13',]

par(mfrow = c(3,2), mar = c(3,3,1,1), oma = c(3,3,0.5,0.5))

for(ispp in 1:(spp + 1)){
  with(subset(x = res_df, spp_scen == ispp),
       plot(tot_mean ~ tot_var, 
            type = 'n', las = 1,
            xlim = range(tot_var), ylim = range(tot_mean) * c(0.9, 1),
            ann = F)) 
  legend('topleft', spp_labels[ispp], bty = 'n', cex = 2)
  
  for(i in c(75, seq(from=50, to=300, by=50))){
    lines(tot_mean ~ tot_var, data = res_df, lwd = 2, 
          subset = (spp_scen == ispp) & (n == i) )
    points(tot_mean ~ tot_var, data = res_df, pch = 16, 
           subset = (spp_scen == ispp) & (n == i))
    with(subset(res_df, (spp_scen == ispp) & (n == i)),
         text(min(tot_var), min(tot_mean), paste('n =', i), pos = 1)
    )
  }
}

mtext(side = 1, 'Variance Criterion', outer = T, line = 1, font = 2)
mtext(side = 2, 'Total Mean Criterion', outer = T, line = 1, font = 2)

###########################
## Show specific solutions for a given effort level
## Solution with highest mean
## Solution with highest variance
## Compromise solution 
###########################
for(ispp in 1:(spp + 1)){
  par(mfrow = c(4,4), oma = c(0,0,3,0))
  
  for(i in c(100, 150, 200, 250)){
    temp_df = subset(res_df, n == i & spp_scen == ispp)
    temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
                       ncol = nrow(stations) )
    
    max_mean_idx = which.max(temp_df$tot_mean)
    max_var_idx = which.max(temp_df$tot_var)
    
    range_mean = diff(range(temp_df$tot_mean))
    range_var = diff(range(temp_df$tot_var))
    
    dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
    dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
    
    comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
    
    sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
    
    par(mar = c(5,4,1,1))
    plot(tot_mean ~ tot_var, data = temp_df, 
         type = 'b', las = 1, pch = 16, lwd = 2, cex = 2, col = 'grey',
         xlim = range(temp_df$tot_var), ylim = range(temp_df$tot_mean),
         xlab = 'Total Variance', ylab = "Total Mean")
    legend('bottomleft', paste('n =', i), bty = 'n', cex = 2, text.font = 2)
    points(tot_mean ~ tot_var, data = temp_df[sol_idx,], pch = 16, cex = 3, 
           col = c('black', 'darkblue', 'red'))
    
    par(mar = rep(0.5, 4))
    for(s in 1:3){
      plot(lat ~ long, data = stations, ann = F, axes = F, asp = 1); box()
      points(lat ~ long, data = stations[temp_mat[sol_idx[s],] == 1,], 
             pch = 16, cex = 2, col = c('black', 'darkblue', 'red')[s])
    }
  }
  
  mtext(side = 3, spp_labels[ispp], outer = T)
  
}
