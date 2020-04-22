# shape_dir = 'C:\\Users\\Zack Oyafuso\\Desktop\\survey_grids\\'
# goa = readOGR(paste0(shape_dir, 'goagrid2019_landuntrawlsndmn.shp'))

########################
## Build and Run Multispecies VAST
########################

rm(list = ls())

library(VAST); library(sp); library(rgdal)

github_dir = 'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/VAST_output8/'

## Import Settings
load(paste0(github_dir, 'Model_Settings.RData'))
load(paste0(github_dir, 'Spatial_Settings.RData'))

# Build and run model

## Build model: To estimate parameters, we first build a list of data-inputs used for parameter estimation.  `make_data` has some simple checks for buggy inputs, but also please read the help file `?make_data`.  
#We then build the TMB object.

TmbData = make_data("Version"=Version, 
                    "FieldConfig"=FieldConfig, 
                    "OverdispersionConfig"=OverdispersionConfig, 
                    "RhoConfig"=RhoConfig, 
                    "ObsModel_ez" = c(PosDist = 2, Link = 0), 
                    "c_i"=as.numeric(Data_Geostat[,'spp'])-1, 
                    "b_i"=Data_Geostat[,'Catch_KG'], 
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
                    "s_i"=Data_Geostat[,'knot_i']-1, 
                    "t_i"=Data_Geostat[,'Year'], 
                    "spatial_list"=Spatial_List, 
                    "Options"=Options ,
                    formula = "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                    covariate_data = cbind(Data_Geostat[,c('Lat', 'Lon', 
                                                           'LOG_DEPTH',
                                                           'LOG_DEPTH2',
                                                           'Catch_KG')], 
                                           Year = NA)
)

str(TmbData$X_gtp)


#Add "true" and not interpolated covariate data
load(paste0(dirname(github_dir), '/Extrapolation_depths.RData'))
X_gtp = array(dim = c(TmbData$n_g,TmbData$n_t, 1) )
for(i in 1:TmbData$n_t) {
  X_gtp[,i,] = as.matrix(Extrapolation_depths[,c('depth')])
}

par(mar = c(4,4,1,1))
plot( x = X_gtp[,1,1], y = TmbData$X_gtp[,1,1], 
      las = 1, pch = 16, cex = 0.25, #log = 'xy',
      xlab = 'Bathymetry from MarMap', ylab = 'Interpolated Survey Bathymetry')
abline(a = 0, b = 1, col = 'red', lwd = 3)

# cor(log10(X_gtp[,1,1]), log10(TmbData$X_gtp[,1,1]))

rect(xleft = max(Data_Geostat$LOG_DEPTH), 
     xright = max(Data_Geostat$LOG_DEPTH)*3,
     ybottom = 0, ytop = 1000, col = hsv(1,1,1,.25))
rect(xleft = -50, 
     xright = min(Data_Geostat$LOG_DEPTH)*3,
     ybottom = 0, ytop = 1000, col = hsv(1,1,1,.25))
