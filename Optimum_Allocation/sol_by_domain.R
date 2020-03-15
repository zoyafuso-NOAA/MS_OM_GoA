load('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Extrapolation_depths.RData')

goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                             data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                          X1 = solution$framenew$STRATO,
                                          domain = df$Domain) )

{tiff('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/solution_map.tiff', res = 200, width = 190, height = 200, units = 'mm',
      compression = 'lzw')
 par(mfrow = c(3,2), oma = rep(1,4), family = 'serif' )
 
 #Plot by management region
 for(idom in 1:ndom){
  goa_ras = raster(subset(goa, domain == levels(goa$domain)[idom] ), resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'X1')
  
  colorRamp = colorRampPalette(
   
   list(c('yellow', 'red', 'brown'),
        c('grey', 'yellow', 'gold'),
        c('lawngreen', 'green', 'darkgreen'),
        c('grey', 'blue', 'darkblue'),
        c('grey', 'purple', 'darkorchid4'))[[idom]]
  )(10)
  
  par(mar = c(0,0,0,0))
  plot(goa_ras, 
       col = colorRamp[c(10,7,3)], 
       axes = F, legend = F )
  
  legend_label = c()
  temp_df = subset(strataStructure, Domain == idom)
  
  for(istratum in 1:nrow(temp_df)){
   legend_label = c( legend_label, paste0( 'Str ', istratum, ': ', temp_df$Population[istratum],
                                           ' units--', temp_df$Allocation[istratum], ' allocated'))
  }
  
  xrange = diff(par()$usr[1:2])
  yrange = diff(par()$usr[3:4])
  legend(x = par()$usr[1]+xrange*c(0.1,0.1,0.45,0.25,0.45)[idom],
         y = par()$usr[3]+yrange*c(0.975,0.975,0.275,0.98,0.975)[idom],
         legend = legend_label, bty = 'n', pt.cex = 3,
         title = paste0('Region ', idom),
         col = colorRamp[c(10,7,3)], pch = 15, cex = 1.)
 }
 
 par( mar = c(0,0,0,0) )
 plot(N_km ~ E_km, data = Extrapolation_depths, pch = '.', cex = 2,
      col = brewer.pal(n = 10, name = 'Spectral')[c(1,5,7,9,10)][frame$domainvalue],
      axes = F, asp=1)
 text(x = tapply(Extrapolation_depths$E_km, frame$domainvalue, mean),
      y = tapply(Extrapolation_depths$N_km, frame$domainvalue, mean),
      font = 2, cex = 2)
 
 xrange = diff(par()$usr[1:2])
 yrange = diff(par()$usr[3:4])
 text(x = par()$usr[1]+xrange*0.35,
      y = par()$usr[3]+yrange*0.85,
      paste0('Optimal Sample Size: ', sum(strataStructure$Allocation), '\n',
             'CV constraint: ', unique(cv$CV1)*100, '%'),
      cex = 2)
 
 #expected_CV(solution$aggr_strata)
 dev.off()
}