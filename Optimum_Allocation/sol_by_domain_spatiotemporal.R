library(raster); library(RColorBrewer); library(SamplingStrata)

which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
result_wd = c('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/', 
              'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')[which_machine]

setwd(result_wd)
load('../Extrapolation_depths.RData')
VAST_model = '6g'
load(paste0('model_', VAST_model, '/optimization_spatiotemporal.RData'))

res_df = as.data.frame(res_df)

strata_list = strata_list[1:2]
sample_sizes = sapply(strata_list, FUN = function(x) sum(x$Allocation))

winner = which.min(sample_sizes)

ns = 15

{tiff(paste0('model_', VAST_model, '/solution_map_spatiotemporal.tiff'), 
      res = 200, width = 8, height = 6, units = 'in', compression = 'lzw')
  
  #Plot by management region
  goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                               data = data.frame(X1=res_df[,paste0('V',winner+2)]) )
  
  goa_ras = raster(goa, resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'X1')
  
  temp_df = strata_list[[winner]]
  nstrata = nrow(temp_df)
  
  par(mar = c(0,0,0,0))
  image(goa_ras, asp = 1, axes = F, col = brewer.pal(n = 12, 'Paired'), ann = F )
  
  legend_label = c()
  
  
  for(istratum in 1:nrow(temp_df)){
    legend_label = c( legend_label, 
                      paste0(temp_df$Population[istratum],
                             ' units--', temp_df$Allocation[istratum], 
                             ' allocated (', temp_df$Lower_X1[istratum], '-',
                             temp_df$Upper_X1[istratum], ' m)'))
  }
  
  xrange = diff(par()$usr[1:2])
  yrange = diff(par()$usr[3:4])
  legend(x = par()$usr[1]+xrange*c(0.1,0.1,0.45,0.25,0.45)[idom],
         y = par()$usr[3]+yrange*c(0.975,0.975,0.275,0.98,0.975)[idom],
         legend = legend_label, bty = 'n', pt.cex = 3,
         title = paste0('Region ', idom),
         col = colorRamp[c(9,7,5,3,1)[1:nstrata]], pch = 15, cex = 0.75)
  
  
  par( mar = c(0,0,0,0) )
  plot(N_km ~ E_km, data = Extrapolation_depths, pch = '.', cex = 2,
       col = brewer.pal(n = 10, name = 'Spectral')[c(1,5,7,9,10)][goa$domainvalue],
       axes = F, asp=1)
  text(x = tapply(Extrapolation_depths$E_km, goa$domainvalue, mean),
       y = tapply(Extrapolation_depths$N_km, goa$domainvalue, mean),
       font = 2, cex = 2)
  
  xrange = diff(par()$usr[1:2])
  yrange = diff(par()$usr[3:4])
  text(x = par()$usr[1]+xrange*0.35,
       y = par()$usr[3]+yrange*0.85,
       paste0('Optimal Sample Size: ', 
              sum(strata_list[[winner]]$Allocation), '\n',
              'CV constraint: ', unique(0.3)*100, '%'),
       cex = 2)
  
  dev.off()
}

temp_df
temp = subset(frame, select = c(id, X1, X2))

depth_cuts = round(sort(unique(c(temp_df$Lower_X1, temp_df$Upper_X1))))
lon_cuts = round(sort(unique(c(temp_df$Lower_X2, temp_df$Upper_X2))))

solution_plane = matrix(nrow = length(depth_cuts)-1, ncol = length(lon_cuts)-1,
                        dimnames = 
                          list(paste0(depth_cuts[1:(length(depth_cuts)-1)], 
                                      '-', 
                                      depth_cuts[2:length(depth_cuts)]),
                               paste0(lon_cuts[1:(length(lon_cuts)-1)], 
                                      '-', 
                                      lon_cuts[2:length(lon_cuts)])))

for(i in 1:(length(lon_cuts)-1) ){
  for(j in 1:(length(depth_cuts)-1) ){
    solution_plane[j,i] = any(temp$X2>lon_cuts[i] & temp$X2<lon_cuts[i+1] & temp$X1>depth_cuts[j] & temp$X1<depth_cuts[j+1]) 
  }
}

solution_plane = solution_plane[!apply(solution_plane, 
                                       MARGIN = 1, 
                                       FUN = function(x) all(x == F))
                                , ]

solution_plane = solution_plane[,!apply(solution_plane, 
                                       MARGIN = 2, 
                                       FUN = function(x) all(x == F)) ]
