###############################################################################
## Project:       Extract Density Values from VAST output
## Author:        Zack Oyafuso
## Description:   For a given VAST model, extract predicted density
###############################################################################

#Import some spatial libaries
library(sp)
library(raster)

#Define working directory of the VAST output files
VAST_dir = 'C:/Users/Zack Oyafuso/Google Drive/GOA_VAST_Runs/VAST_output10a/'
github_dir = 'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/'

#Load VAST output
load(paste0(VAST_dir, 'fit.RData'))

#Load CPUE dataset
load( paste0(dirname(VAST_dir), '/Spatial_Settings_CrVa.RData') )

#Load Spatial Information
load(paste0(github_dir, 'data/Extrapolation_Depths.RData'))

###########################################
## Predicted Density (fit$Report$D_gcy) is indexed by: 
## Extrapolation grid cell (g): 23339
## Species (c): 15 total species (sci_names)
## year (y): 24 years (Year_Set1996-2019), 
##           but only 11 years with data (Years2Include)
###########################################
sci_names = sort(unique(Data_Geostat$spp))

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

str(fit$Report$D_gcy[,,Years2Include])

#An easy way to visualize density distribution
#Sebastes polyspinis in 1996 as an example
plot_this_density = fit$Report$D_gcy[,13,24]

goa = SpatialPointsDataFrame(
  coords = Extrapolation_depths[,c('E_km', 'N_km')],
  data = data.frame(density=plot_this_density) )
goa_ras = raster(goa, resolution = 5)
goa_ras =rasterize(x = goa, y = goa_ras, field = 'density')
plot(goa_ras, col = rev(terrain.colors(1000)), axes = F)
