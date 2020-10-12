###############################################################################
## Project:       VAST diagnostics plots
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate various diagnostic plots
##                Data and knot locations
##                Encounter probability
##                QQ plots
##                Pearson Residuals 
##                Mean density andyear-specific density, 
##                Spatial (omega) and spatiotemporal (epsilon) effects, 
##                Index of abundance
##                Covariance
##                Factors loadings
##                Anisotropy
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries
##################################################
library(VAST)
library(raster)
library(sp)
library(rgdal)
library(RColorBrewer)
library(plotrix)
library(rnaturalearth)
library(rgeos)

##################################################
#### Set up directories     
##################################################
depth_in_model <- c(T, F)[1]
which_spp <- c("Atheresthes stomias", "Gadus chalcogrammus", 
               "Gadus macrocephalus", "Glyptocephalus zachirus",
               "Hippoglossoides elassodon", "Hippoglossus stenolepis",
               "Lepidopsetta bilineata", "Lepidopsetta polyxystra",
               "Limanda aspera", "Microstomus pacificus", "Sebastes alutus",
               "Sebastes B_R", "Sebastes polyspinis", "Sebastes variabilis", 
               "Sebastolobus alascanus" ) [14]
# ns <- length(sci_names)

VAST_dir <- paste0("G:/Oyafuso/VAST_Runs_EFH/Single_Species/", which_spp,
                   ifelse(depth_in_model, "_depth", ""), "/")
fun_dir <- paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/diagnostics/")
diag_dir <-  paste0(VAST_dir, "diagnostics/")

if (! dir.exists(diag_dir)) dir.create(diag_dir)

##################################################
####  Load VAST fit and import customized plot functions 
##################################################
load(paste0(VAST_dir, "/fit.RData"))
# load(paste0(VAST_dir, "/Spatial_Settings_CrVa.RData") )

# for (ifile in c("plot_factors.R",
#                 "plot_residuals.R",
#                 "summarize_covariance.R")) {
#   source( paste0(fun_dir, ifile) )
#   rm(ifile)
# }

##################################################
#### Extract objects from fitted object   
#### Set up constants
##################################################
Data_Geostat <- fit$data_frame
names(Data_Geostat)[c(1:2, 5:7)] <- c("Lat", "Lon", "Catch_KG", "Year", "spp")

Spatial_List <- fit$spatial_list
Extrapolation_List <- fit$extrapolation_list
Opt <- fit$parameter_estimates
Obj <- fit$tmb_list$Obj
Report<- fit$Report
TmbData <- fit$data_list

Year_Set <- seq(min(Data_Geostat[,"Year"]), max(Data_Geostat[, "Year"]))
Years2Include <- which(Year_Set %in% sort(unique(Data_Geostat[, "Year"])))

MapDetails_List <- make_map_info("Region" = "Gulf_of_Alaska", 
                                 "spatial_list" = Spatial_List, 
                                 "Extrapolation_List" = Extrapolation_List )


xrange <- range(Extrapolation_List$Data_Extrap[, "E_km"])
yrange <- range(Extrapolation_List$Data_Extrap[, "N_km"])
xrange_diff <- diff(xrange)
yrange_diff <- diff(yrange)

##################################################
####   Import Land objects, clip to save object space
##################################################
AK = readOGR(paste0(dirname(fun_dir), "/data/shapefiles/AKland.shp")) 
AK = sp::spTransform(AK, 
                     CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                  " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                  " +ellps=GRS80 +datum=NAD83 +units=km",
                                  " +no_defs"))

CA = readOGR(paste0(dirname(fun_dir), "/data/shapefiles/canada_dcw.shp")) 
CA = sp::spTransform(CA, 
                     CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                  " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                  " +ellps=GRS80 +datum=NAD83 +units=km",
                                  " +no_defs"))

AK_bbox = matrix(c(-886, 2163, 
                   5600, 7000), 
                 byrow = T, ncol = 2)
AK = gIntersection(AK, 
                   as(extent(as.vector(t(AK_bbox))), "SpatialPolygons"), 
                   byid = TRUE)
AK = rgeos::gSimplify(spgeom = AK, tol = 3)

CA = subset(CA, POPYADMIN %in% c("BRITISH COLUMBIA",
                                 "YUKON TERRITORY"))
CA = rgeos::gSimplify(spgeom = CA, tol = 3)

############################################
## If continuting from where you left on...
############################################
# load(paste0(diag_dir, "/diagnostics.RData"))

############################################
## Plot data
## It is always good practice to conduct exploratory analysis of data.  Here, I
## visualize the spatial distribution of data.  Spatio-temporal models involve 
## the assumption that the probability of sampling a given location is 
## statistically independent of the probability distribution for the response 
## at that location.  So if sampling "follows" changes in density, then the 
## model is probably not appropriate!
############################################

# plot_data(Extrapolation_List = Extrapolation_List,
#           Spatial_List = Spatial_List,
#           Data_Geostat = Data_Geostat,
#           Year_Set = Years2Include,
#           PlotDir = diag_dir)

############################################
## Convergence
## Here I print the diagnostics generated during parameter estimation, and I 
## confirm that (1) no parameter is hitting an upper or lower bound and (2) the
## final gradient for each fixed-effect is close to zero. For explanation of 
## parameters, please see `?make_data`.
############################################

# pander::pandoc.table( Opt$diagnostics[,c("Param","Lower","MLE",
#                                          "Upper","final_gradient")] ) 
# 
# max(abs(Opt$diagnostics$final_gradient))

############################################
## Diagnostics for encounter-probability component
## Next, we check whether observed encounter frequencies for either low or high 
## probability samples are within the 95% predictive interval for predicted 
## encounter probability
############################################
Enc_prob = plot_encounter_diagnostic( Report = Report, 
                                      Data_Geostat = Data_Geostat, 
                                      DirName = diag_dir)

############################################
## Diagnostics for positive-catch-rate component
## We can visualize fit to residuals of catch-rates given encounters using a
## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
############################################

# Q = plot_quantile_diagnostic( TmbData = TmbData,
#                               Report = Report,
#                               FileName_PP = "Posterior_Predictive",
#                               FileName_Phist = "Posterior_Predictive-Histogram",
#                               FileName_QQ = "Q-Q_plot",
#                               FileName_Qhist = "Q-Q_hist",
#                               DateFile = diag_dir )

############################################
## Save output
############################################
# save(list = "Q",
#      file = paste0(diag_dir, "diagnostics.RData"))

############################################
## Diagnostics for positive-catch-rate component
## We can visualize fit to residuals of catch-rates given encounters using a
## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
############################################

# {png(filename = paste0(diag_dir, "QQplot.png"), 
#      width = 6, 
#      height = 8, 
#      units = "in", 
#      res = 200)
#   
#   par(mfrow = c(5, 3), 
#       mar = c(0, 0, 0, 0), 
#       oma = c(5, 5, 1, 1))
#   
#   for (ispp in 1:ns) {
#     Q_temp = na.omit(Q[[ispp]]$Q)
#     Order = order(Q_temp)
#     plot(x = seq(0, 1, length = length(Order)), y = Q_temp[Order], 
#          ylim = c(0,1.0), xlim = c(0,1.0), 
#          type = "l", lwd = 3, las = 1, axes = F, ann = F)
#     abline(a = 0, b = 1)
#     box()
#     text(0.45, 0.950, sci_names[ispp], font = 3)
#     
#     if(ispp%%6 == 1) axis(side = 2, las = 1)
#     if(ispp %in% c(ns-2, ns) ) axis(side = 1)
#   }
#   rm(Q_temp, ispp, Order)
#   mtext(side = 1, "Theoretical Quantile", outer = T, line = 3)
#   mtext(side = 2, "Sample Quantile", outer = T, line = 3)
#   dev.off()}

############################################
## We then plot Pearson residuals.  If there are visible patterns (areas with 
## consistently positive or negative residuals accross or within years) then 
## this is an indication of the model "overshrinking" results towards the 
## intercept, and model results should then be treated with caution.  
############################################
# if(!dir.exists(paste0(diag_dir,"Pearson_Residuals/"))) 
#   dir.create(paste0(diag_dir,"Pearson_Residuals/"))
# 
# PResid = plot_residuals(Lat_i=Data_Geostat[,"Lat"], 
#                         Lon_i=Data_Geostat[,"Lon"], 
#                         TmbData=TmbData, 
#                         Report=Report, 
#                         Q=Q, 
#                         spatial_list = Spatial_List,
#                         extrapolation_list = Extrapolation_List,
#                         working_dir=paste0(diag_dir,"Pearson_Residuals/"), 
#                         Year_Set=Year_Set, 
#                         Years2Include=Years2Include, 
#                         mar=c(0,0,2,0), 
#                         oma=c(3.5,3.5,0,0), 
#                         cex=1.8)
# 
# save(list = c("Q", "PResid"),
#      file = paste0(diag_dir, "diagnostics.RData"))

############################################
## Plot mean Density for each Species
############################################
# {
#   #Density colors
#   colors = c("white", brewer.pal(n = 9, "Oranges"), "black")
#   
#   png(filename = paste0(PP_dir, "mean_annual_density_", modelno, ".png"), 
#       units = "in", height = 6, width = 9, res = 500)
#   
#   ## Set up plot panel layots
#   par(mar = c(0,0,0,0))
#   layout(mat = matrix(c(1,2,3,16,
#                         4,5,6,16,
#                         7,8,9,16,
#                         10,11,12,16,
#                         13,14,15,16), byrow = T, nrow = 5), 
#          widths = c(rep(1,3),0.4) )
#   
#   for(ispp in 1:ns){
#     
#     ##Empty plot
#     plot(1, type = "n", axes = F, ann = F, xlim = xrange, ylim = yrange )
#     
#     #Extract mean density across years for each sepcies
#     mean_dens_years  = rowMeans(Report$D_gcy[,ispp,Years2Include])
#     goa = SpatialPointsDataFrame(
#       coords = Extrapolation_List$Data_Extrap[,c("E_km", "N_km")], 
#       data = data.frame(var = mean_dens_years))
#     goa_ras = raster(goa, resolution = 5)
#     goa_ras = rasterize(x = goa, y = goa_ras, field = "var")
#     
#     #Tabularize density values by deciles
#     val_cuts = c(0,1, quantile(mean_dens_years, probs = seq(0,1,0.1)))
#     values(goa_ras) = cut(x = values(goa_ras),
#                           breaks = val_cuts)
#     
#     #Plot density
#     image(goa_ras, asp = 1, axes = F, ann = F, add = T, col = colors, asp = 1)
#     
#     #Species Label
#     text(x = xrange[1] + xrange_diff*0.71,
#          y = yrange[1] + yrange_diff*0.25,
#          gsub(sci_names[ispp], pattern = " ", replacement = "\n"), 
#          cex = 1.25, font = 3)
#     
#     #Add land
#     plot(AK, add = T, col = "tan", border = F)
#     plot(CA, add = T, col = "tan", border = F)
#     
#     box()
#   }
#   
#   #Plot legend
#   plot(1, type = "n", axes = F, ann = F, xlim = c(0,1), ylim = c(0,10))
#   plotrix::color.legend(xl = 0.05,
#                         xr = 0.35,
#                         yb = 1,
#                         yt = 9,
#                         legend = paste0(seq(0,100,10), " %" ),
#                         rect.col = colors,
#                         gradient = "y", align = "rb")
#   dev.off()
# }


############################################
## Plot Density across years for each Species
############################################
{
  pdf(paste0(diag_dir, "density.pdf"), 
      height = 5, width = 7)
  
  #Plot layout
  par(mar = c(0,0,0,0), oma = rep(0.5,4), mfrow = c(4,3))
  
  for(iyear in Years2Include){
    
    #Extract density values for a species in a year, 
    vals  = Report$D_gct[,1,iyear]
    val_cuts = c(0,quantile(vals[vals > 1], probs = seq(0,1,length=9) ))
    
    #plot density
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_List$Data_Extrap[,c("E_km", "N_km")], 
      data = data.frame(density = vals) )
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = "density")
    
    values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
    
    colors = c("white", brewer.pal(n = 7, name = "Oranges"), "black")
    image(goa_ras, asp = 1, axes = F, ann = F, add = F, 
          col = colors)
    
    #Year label
    text(x = goa_ras@extent[1] + 0.7*diff(goa_ras@extent[1:2]),
         y = goa_ras@extent[3]+ 0.7*diff(goa_ras@extent[3:4]),
         Year_Set[iyear], cex = 1)
    
    #Add value legend
    val_cuts = round(val_cuts[-1])
    legend("bottom", fill = colors, bty = "n",
           ncol = 3, cex = 0.65,
           legend = c("<1", paste0("1-", val_cuts[2]), 
                      paste0(val_cuts[2:(length(val_cuts)-1)], "-",
                             val_cuts[3:length(val_cuts)])) )
    
    #Add land
    plot(AK, add = T, col = "tan", border = F)
    plot(CA, add = T, col = "tan", border = F)
    
    box()  
  }
  plot(1, type = "n", axes = F, ann = F)
  caption = paste0("Predicted density (kg/km2) for\n", which_spp,
                   "\nacross the Gulf of Alaska\nfor each observed year.")
  text(1,1, caption, cex = 1.2, font = 3)
  
  dev.off()
}

############################################
## Plot Omega of the 1st and second components
############################################
{
  png(filename = paste0(diag_dir, "omega.png"),
      units = "in", height = 4, width = 6, res = 500)
  
  #Plot layout
  par(mar = c(0,0,0,0), oma = c(0,0,2,0))
  layout(mat = matrix(c(1,2), byrow = T, nrow = 1), 
         widths = c(1, 0.25) )
  
  
  offset = 0
  
  #Empty plot
  plot(1, type = "n", axes = F,
       xlim = xrange, ylim = yrange + c(-0.5*yrange_diff, 0))
  for(omegatype in 1:2){ #Two types for the 0/1 and pos components
    
    #Extract spatial component
    scaled_var = list(Report$Omega1_gc,
                      Report$Omega2_gc)[[omegatype]]
    
    #Scale to standard normal
    scaled_var = (scaled_var - mean(scaled_var)) / sd(scaled_var)
    
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_List$Data_Extrap[,c("E_km", "N_km")], 
      data = data.frame(var = scaled_var) )
    
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = "var")
    
    goa_ras = raster::shift(goa_ras, dy = -offset*yrange*0.075)
    
    #Plot spatial effect
    colors = rev(brewer.pal(n = 11, name = "Spectral"))
    image(goa_ras, add = T, col = colors, asp = 1)
    
    temp_yrange = extent(goa_ras)[3:4]
    if(omegatype == 2) {
      #Plot Species label
      text(x = xrange[1] + xrange_diff*0.725,
           y = temp_yrange[1] + yrange_diff*0.3,
           gsub(x = which_spp, pattern = " ", replacement = "\n"), 
           font = 3, cex = 1.25)
    }
    
    offset=offset + 1
    box()
  }
  
  
  #Plot legend
  plot(1, type = "n", axes = F, ann = F, xlim = c(0,1), ylim = c(0,10))
  plotrix::color.legend(xl = 0.1,
                        xr = 0.5,
                        yb = 1,
                        yt = 9,
                        legend = -3:3,
                        rect.col = colorRampPalette(colors)(1000) ,
                        gradient = "y", align = "rb")
  
  #Title
  mtext(side = 3, outer = T, 
        "Spatial Effect in Occurrence (Top) and Positive Response (Bottom)", 
        line = 0, font = 2)
  dev.off()
  rm(offset, omegatype, scaled_var, temp_yrange)
}

############################################
## Plot Epsilon for the first, middle and last year
## for each of the 1st and second components
############################################
png(paste0(diag_dir, "Epsilon.png"),
    units = "in", 
    height = 5, 
    width = 6, 
    res = 500)

#Plot layout
par(mar = c(0,0,0,0), oma = c(0.5, 0.5, 3, 0.5), mfrow = c(4,3))
for(iyear in Years2Include){
  
  #Empty plot
  plot(1, type = "n",  axes = F,
       xlim = xrange, ylim = yrange + c(-1*yrange_diff, 0))
  
  #Year label
  legend("topleft", legend = Year_Set[iyear], bty = "n")
  box()
  offset = 0
  
  for(epstype in 1:2){ #Two types for the 0/1 and positive components
    
    #Extract spatiotemporal component
    scaled_var = list(Report$Epsilon1_gct[,,iyear],
                      Report$Epsilon2_gct[,,iyear])[[epstype]]
    
    #Scale to standard normal
    scaled_var = (scaled_var - mean(scaled_var)) / sd(scaled_var)
    
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_List$Data_Extrap[,c("E_km", "N_km")],
      data = data.frame(var = scaled_var) )
    
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = "var")
    
    goa_ras = raster::shift(goa_ras, dy = -offset*yrange*0.12)
    
    #Plot spatiotemporal effect
    colors = rev(brewer.pal(n = 11, name = "Spectral"))
    image(goa_ras, add = T, col = colors, asp = 1)
    
    temp_yrange = extent(goa_ras)[3:4]
    offset=offset + 1
  }
}

#Plot legend
plot(1, type = "n", axes = F, ann = F, xlim = c(0,1), ylim = c(0,10))
plotrix::color.legend(xl = 0.05, xr = 0.95, yb = 3.5, yt = 5.5, align = "rb",
                      legend = -3:3, gradient = "x", cex = 0.75, 
                      rect.col = colorRampPalette(colors)(1000))

#Species label
mtext(side = 3, which_spp, line = -3, font = 3, cex = 1)

#Title
mtext(side = 3, outer = T, 
      paste0("Spatiotemporal Effect in Occurrence (Top)", 
             " and Positive Response (Bottom)"), line = 1)


dev.off()


############################################
## Index of abundance
############################################

Index = plot_biomass_index( DirName=diag_dir, 
                            TmbData=TmbData, 
                            Sdreport=Opt[["SD"]], 
                            Year_Set=Year_Set, 
                            Years2Include=Years2Include, 
                            strata_names = strata.limits <- data.frame(
                              "STRATA" = c("All_areas"),#, "west_of_140W"),
                              "west_border" = c(-Inf),#, -Inf),
                              "east_border" = c(Inf)#, -140)
                            )[,1], 
                            use_biascorr=TRUE, 
                            category_names=levels(Data_Geostat[,"spp"]) )

############################################
## Plot spatial and spatio-temporal covariance
## We can visualize the spatial and spatio-temporal covariance among species 
## in encounter probability and positive catch rates (depending upon what is
### turned on via `FieldConfig`):
############################################
# Cov_List = summarize_covariance( Report=Report, 
#                                  ParHat=Obj$env$parList(), 
#                                  Data=TmbData, 
#                                  SD=Opt$SD, 
#                                  plot_cor=FALSE, 
#                                  category_names=levels(Data_Geostat[,"spp"]), 
#                                  plotdir=diag_dir, 
#                                  mgp=c(2,0.5,0), 
#                                  tck=-0.02, 
#                                  oma=c(0,10,2,2) )

############################################
## Plot factors
## Finally, we can inspect the factor-decomposition for community-level 
## patterns.  This generates many plots, only some of which are included in 
## this tutorial document.
############################################
# if(!dir.exists(paste0(diag_dir,"Factors/"))) 
#   dir.create(paste0(diag_dir,"Factors/"))
# 
# Plot_factors( Report=Report, 
#               ParHat=Obj$env$parList(), 
#               Data=TmbData, 
#               SD=Opt$SD, 
#               mapdetails_list=MapDetails_List, 
#               Year_Set=Year_Set, 
#               category_names=levels(Data_Geostat$spp), 
#               plotdir=paste0(diag_dir,"Factors/") )

############################################
## Direction of "geometric anisotropy"
## We can visualize which direction has faster or slower decorrelation 
## (termed "geometric anisotropy")
############################################
# plot_anisotropy( FileName=paste0(diag_dir,"Aniso.png"), 
#                  Report=Report, 
#                  TmbData=TmbData )

############################################
## Save output
############################################
# save(list = c("Q", "PResid", "Cov_List"),
#      file = paste0(diag_dir, "diagnostics.RData"))

