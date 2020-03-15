library(raster); library(RColorBrewer); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')
load('../Extrapolation_depths.RData')

VAST_model = '6d'
load(paste0('model_', VAST_model, '/optimization.RData'))
res_df = as.data.frame(res_df)


sample_sizes = sapply(strata_list, FUN = function(x) sum(x$Allocation))

winner = which.min(sample_sizes)

ns = 15
goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                             data = data.frame(
                               id = res_df$id,
                               domainvalue = res_df$domainvalue,
                               X1 = res_df[,paste0('V',winner)] ) )

# goa = SpatialPoints(coords = Extrapolation_depths[,c('E_km', 'N_km')] )
goa_ras = raster(goa, resolution = 5)
goa_ras =rasterize(x = goa, y = goa_ras, field = 'X1')

par(mar = c(0,0,0,0))
nstrata = nrow(strata_list[[winner]])

image(goa_ras, col = 'white', asp = 1, axes = F, ann = F)

for(idom in 1:5){
  sub_goa = subset(goa, domainvalue == idom)
  sub_goa_ras = raster(sub_goa, resolution = 5)
  sub_goa_ras = rasterize(x = subset(goa, domainvalue == idom), y = sub_goa_ras, field = 'X1')
  
  nstrata = length(unique(subset(goa, domainvalue == idom)$X1))
  ncolors = brewer.pal(n = nstrata, name = 'Spectral')
  
  plot(sub_goa_ras, add = T, col = ncolors, legend = F)
}

sum(strata_list[[winner]]$Allocation)

