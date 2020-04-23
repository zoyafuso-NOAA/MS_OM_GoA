##################################
## Diagnostic and Output plots
##################################

rm(list = ls())

library(VAST)
library(raster)
library(sp)
library(RColorBrewer)
library(plotrix)
# library(FishStatsUtils)

modelno_main = '7'
modelno = '7a'

VAST_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output", 
                  modelno_main, "/VAST_output", modelno, "/")
diag_dir = paste0("C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/diagnostics/",
                  "VAST_model", modelno_main, "/VAST_output", modelno, "/")
fun_dir = paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/diagnostics/')

#################################
## Load Data and modified functions
#################################

for(ifile in c("plot_factors.R", "plot_maps_density.R",
               "plot_variable_density.R", "summarize_covariance.R",
               'plot_variable.R')){
  source( paste0(fun_dir, ifile) )
  rm(ifile)
}

# load( paste0(dirname(VAST_dir), '/Spatial_Settings_CrVa.RData') )
load(paste0(diag_dir, 'diagnostics.RData'))

###################################
## Plot Settings
###################################
#Dummy variable just to extract ranges
depth_scale = raster(plot_variable(Y_gt = TmbData$X_gtp[,1,1],
                                   map_list = MapDetails_List,
                                   working_dir = diag_dir,
                                   file_name = 'depth',
                                   oma = c(0,0,0,0), mar = c(1,1,1,1) ,
                                   projargs = "+proj=utm +zone=5")[[1]])

xrange = range(extent(depth_scale)[1:2])
yrange = range(extent(depth_scale)[3:4])
xrange_diff = diff(xrange)
yrange_diff = diff(yrange)

##################################
## Extract objects from fitted object
##################################

######################
## Plot data
## It is always good practice to conduct exploratory analysis of data.  Here, I 
## visualize the spatial distribution of data.  Spatio-temporal models involve 
## the assumption that the probability of sampling a given location is 
## statistically independent of the probability distribution for the response 
## at that location.  So if sampling "follows" changes in density, then the 
## model is probably not appropriate!
############################################

# plot_data(Extrapolation_List=Extrapolation_List, 
#           Spatial_List=Spatial_List, 
#           Data_Geostat=Data_Geostat, 
#           Year_Set = Years2Include,
#           PlotDir=diag_dir)

############################################
## Convergence
## Here I print the diagnostics generated during parameter estimation, and I 
## confirm that (1) no parameter is hitting an upper or lower bound and (2) the
## final gradient for each fixed-effect is close to zero. For explanation of 
## parameters, please see `?make_data`.
############################################

pander::pandoc.table( Opt$diagnostics[,c('Param','Lower','MLE',
                                         'Upper','final_gradient')] ) 

max(abs(Opt$diagnostics$final_gradient))

############################################
## Diagnostics for encounter-probability component
## Next, we check whether observed encounter frequencies for either low or high 
## probability samples are within the 95% predictive interval for predicted 
## encounter probability
############################################
Enc_prob = plot_encounter_diagnostic( Report=Report, 
                                      Data_Geostat=Data_Geostat, 
                                      DirName=diag_dir)

############################################
## Diagnostics for positive-catch-rate component
## We can visualize fit to residuals of catch-rates given encounters using a
## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
############################################

{png(filename = paste0(diag_dir, 'QQplot.png'), width = 6, height = 8, 
     units = 'in', res = 200)
  par(mfrow = c(5,3), mar = c(0,0,0,0), oma = c(5,5,1,1))
  for(ispp in 1:ns){
    Q_temp = na.omit(Q[[ispp]]$Q)
    Order = order(Q_temp)
    plot(x = seq(0, 1, length = length(Order)), y = Q_temp[Order], 
         ylim = c(0,1.0), xlim = c(0,1.0), 
         type = "l", lwd = 3, las = 1, axes = F, ann = F)
    abline(a = 0, b = 1)
    box()
    text(0.45, 0.950, sci_names[ispp], font = 3)
    
    if(ispp%%6 == 1) axis(side = 2, las = 1)
    if(ispp %in% c(ns-2, ns) ) axis(side = 1)
  }
  rm(Q_temp, ispp, Order)
  mtext(side = 1, 'Theoretical Quantile', outer = T, line = 3)
  mtext(side = 2, 'Sample Quantile', outer = T, line = 3)
  dev.off()}

############################################
## Plot Covariates
############################################
{
  png(filename = paste0(diag_dir, 'covariates.png'),
      units = 'in', height = 6, width = 6, res = 500)
  
  depth_scale = raster(plot_variable(Y_gt = TmbData$X_gtp[,1,1],
                                     map_list = MapDetails_List,
                                     working_dir = diag_dir,
                                     file_name = 'depth',
                                     oma = c(0,0,0,0), mar = c(1,1,1,1) ,
                                     projargs="+proj=utm +zone=5")$Raster_proj)
  
  depth2_scale = raster(plot_variable(Y_gt = TmbData$X_gtp[,1,2],
                                      map_list = MapDetails_List,
                                      working_dir = diag_dir,
                                      file_name = 'depth',
                                      oma = c(0,0,0,0), mar = c(1,1,1,1) ,
                                      projargs="+proj=utm +zone=5")$Raster_proj)
  
  par(mfrow = c(1,1), mar = c(0,0,0,0))
  plot(1, type = 'n', xlim = xrange, ylim = yrange + c(-3*yrange_diff, 0),
       axes = F, ann = F)
  
  image(depth_scale, asp = 1, col = brewer.pal(n = 10,name = 'Spectral'))
  plotrix::color.legend(xl = xrange[1] + xrange_diff*0.01,
                        xr = xrange[1] + xrange_diff*0.35,
                        yb = yrange[1] + yrange_diff*0.5,
                        yt = yrange[1] + yrange_diff*0.55,
                        legend = -5:3, cex = 0.75,
                        rect.col = brewer.pal(n = 9,name = 'Spectral'),
                        gradient = 'x')#, align = 'rb')
  text(x = xrange[1] + xrange_diff*0.725,
       y = yrange[1] + yrange_diff*0.5,
       'Centered Log-Depth', 
       font = 3, cex = 1.25)
  box(col = 'white')
  
  depth2_scale = raster::shift(depth2_scale, dy = -yrange*0.1)
  image(depth2_scale, asp = 1, add = T, 
        col = brewer.pal(n = 9, name = 'Blues'))
  plotrix::color.legend(xl = xrange[1] + xrange_diff*0.01,
                        xr = xrange[1] + xrange_diff*0.35,
                        yb = yrange[1] - yrange_diff*0.25,
                        yt = yrange[1] - yrange_diff*0.2,
                        legend = seq(0,25,5), cex = 0.75,
                        rect.col = brewer.pal(n = 9, name = 'Blues'),
                        gradient = 'x')#, align = 'rb')
  text(x = xrange[1] + xrange_diff*0.725,
       y = yrange[1] - yrange_diff*0.2,
       'Square of the\nCentered Log-Depth', 
       font = 3, cex = 1.25)
  
  dev.off()
}

############################################
## Plot Density across years for each Species
############################################

{png(paste0(diag_dir, 'density.png'), units = 'in', height = 7, width = 12, res = 500)
  par(mar = c(0,0,1,0), mfrow = c(2,8))
  for(ispp in 1:ns){
    
    vals  = Report$D_gcy[,ispp,]
    val_cuts = c(0,quantile(vals[vals > 1], probs = seq(0,1,0.2) ))
    
    yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
    plot(1, type = 'n', axes = F, ann = F,
         xlim = range(Extrapolation_List$Data_Extrap[,c('E_km')]),
         ylim = c(min(Extrapolation_List$Data_Extrap[,c('N_km')])-5.55*yrange,
                  max(Extrapolation_List$Data_Extrap[,c('N_km')]))
    )
    
    offset = 0
    for(iyear in Years2Include){
      data = as.data.frame(Report$D_gcy[,ispp,iyear])
      names(data) = 'density'
      
      goa = SpatialPointsDataFrame(coords = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], 
                                   data = data )
      goa_ras = raster(goa, resolution = 5)
      goa_ras = rasterize(x = goa, y = goa_ras, field = 'density')
      
      values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
      
      #offset
      goa_ras = raster::shift(x = goa_ras, dy = -yrange/1.85*offset)
      offset = offset + 1
      
      colors = c('white', 'yellow', 'gold', 'orange', 'red', 'brown')
      
      image(goa_ras, asp = 1, axes = F, ann = F, add = T, 
            col = colors)
      
      text(x = goa_ras@extent[1] + 0.7*diff(goa_ras@extent[1:2]),
           y = goa_ras@extent[3]+ 0.7*diff(goa_ras@extent[3:4]),
           Year_Set[iyear], cex = 0.75)
      
    }
    
    mtext(side = 3, gsub(sci_names[ispp], pattern = ' ', replacement = '\n'), 
          line = -1, cex = 0.6, font = 3)
    
    val_cuts = round(val_cuts[-1])
    legend('bottom', fill = colors, bty = 'n',
           ncol = 3, cex = 0.5,
           legend = c('<1', paste0('1-', val_cuts[2]), 
                      paste0(val_cuts[2:(length(val_cuts)-1)], '-',
                             val_cuts[3:length(val_cuts)])) )
  }  
  
  dev.off()}

############################################
## Plot Omega of the 1st and second components
############################################
{
  png(filename = paste0(diag_dir, 'omega.png'),
      units = 'in', height = 6, width = 12, res = 500)
  
  par(mar = c(0,0,0,0), oma = c(0,0,2,0))
  
  layout(mat = matrix(c(1,2,3,4,5, 16,
                        6,7,8,9,10,16,
                        11,12,13,14,15,16), byrow = T, nrow = 3), widths = c(rep(1,5),0.25) )
  
  for(ispp in 1:ns){
    offset = 0
    plot(1, type = 'n', xlim = xrange, ylim = yrange + c(-0.5*yrange_diff, 0),
         axes = F)
    for(omegatype in 1:2){
      scaled_var = list(Report$Omega1_gc[,ispp],
                        Report$Omega2_gc[,ispp])[[omegatype]]
      
      scaled_var = (scaled_var - mean(scaled_var)) / sd(scaled_var)
      
      temp_ras = raster( plot_variable(Y_gt = scaled_var,
                                       map_list = MapDetails_List,
                                       projargs = "+proj=utm +zone=5")[[1]] )
      temp_ras = raster::shift(temp_ras, dy = -offset*yrange*0.075)
      image(temp_ras, add = T,
            col = colorRampPalette(c('blue', 'white', 'red'))(1000))
      
      temp_yrange = extent(temp_ras)[3:4]
      
      if(omegatype == 2) {
        text(x = xrange[1] + xrange_diff*0.725,
             y = temp_yrange[1] + yrange_diff*0.3,
             gsub(x = sci_names[ispp], pattern = ' ', replacement = '\n'), 
             font = 3, cex = 1.25)
      }
      
      offset=offset + 1
      box()
    }
  }
  
  plot(1, type = 'n', axes = F, ann = F, xlim = c(0,1), ylim = c(0,10))
  plotrix::color.legend(xl = 0.1,
                        xr = 0.5,
                        yb = 1,
                        yt = 9,
                        legend = -3:3,
                        rect.col = colorRampPalette(c('blue', 'white', 
                                                      'red'))(1000),
                        gradient = 'y', align = 'rb')
  
  mtext(side = 3, outer = T, 
        "Spatial Effect in Occurrence (Top) and Positive Response (Bottom)", 
        line = 0, font = 2)
  dev.off()
  rm(ispp, offset, omegatype, scaled_var, temp_ras, temp_yrange)
}

############################################
## Plot Epsilon for the first, middle and last year
## for each of the 1st and second components
############################################

{png(filename = paste0(diag_dir, 'epsilon.png'),
     units = 'in', height = 6, width = 12, res = 500)
  layout(mat = matrix(c(1:10,  31,
                        11:20, 31,
                        21:30, 31), byrow = T, nrow = 3), 
         widths = c(rep(1,10),0.25) )
  par(mar = c(1,0,1,0))
  for(ispp in 1:ns){
    for(epstype in 1:2){
      plot(1, type = 'n', xlim = xrange, ylim = yrange + c(-5.5*yrange_diff, 0),
           axes = F)
      offset = 0
      for(iyear in Years2Include){
        scaled_var = list(Report$Epsilon1_gct[,ispp,iyear],
                          Report$Epsilon2_gct[,ispp,iyear])[[epstype]]
        
        scaled_var = (scaled_var - mean(scaled_var)) / sd(scaled_var)
        
        temp_ras = raster( plot_variable(Y_gt = scaled_var,
                                         map_list = MapDetails_List,
                                         projargs = "+proj=utm +zone=5")[[1]] )
        temp_ras = raster::shift(temp_ras, dy = -offset*yrange*0.075)
        image(temp_ras, add = T,
              col = colorRampPalette(c('blue', 'white', 'red'))(1000))
        
        temp_yrange = extent(temp_ras)[3:4]
        
        text(x = xrange[1] + xrange_diff*0.725,
             y = temp_yrange[1] + yrange_diff*0.7,
             Year_Set[iyear], cex = 0.4)
        
        offset=offset + 1
      }
      
      mtext(side = 1, line = 0, cex = 0.4,
            text = c('Spatiotemporal Effect\non Occurrence',
                     'Spatiotemporal Effect\non Positive Response')[epstype])
    }
    text(x = xrange[1], y = yrange[2] + yrange_diff*0.15, sci_names[ispp], 
         font = 3, xpd = NA)
  }
  
  plot(1, type = 'n', axes = F, ann = F, xlim = c(0,1), ylim = c(0,10))
  plotrix::color.legend(xl = 0.05,
                        xr = 0.4,
                        yb = 1,
                        yt = 9,
                        legend = -3:3,
                        rect.col = colorRampPalette(c('blue', 'white', 
                                                      'red'))(1000),
                        gradient = 'y', align = 'rb')
  
  dev.off()}

############################################
## Index of abundance
############################################
if(!dir.exists(paste0(diag_dir,'Index/'))) 
  dir.create(paste0(diag_dir,'Index/'))
Index = plot_biomass_index( DirName=paste0(diag_dir,'Index/'), 
                            TmbData=TmbData, 
                            Sdreport=Opt[["SD"]], 
                            Year_Set=Year_Set, 
                            Years2Include=Years2Include, 
                            strata_names = strata.limits <- data.frame(
                              'STRATA' = c("All_areas"),#, "west_of_140W"),
                              'west_border' = c(-Inf),#, -Inf),
                              'east_border' = c(Inf)#, -140)
                            )[,1], 
                            use_biascorr=TRUE, 
                            category_names=levels(Data_Geostat[,'spp']) )
pander::pandoc.table( Index$Table[,c("Category","Year",
                                     "Estimate_metric_tons","SD_mt")] ) 

######################################
##
######################################

############################################
## Plot factors
## Finally, we can inspect the factor-decomposition for community-level 
## patterns.  This generates many plots, only some of which are included in 
## this tutorial document.
############################################
if(!dir.exists(paste0(diag_dir,'Factors/'))) 
  dir.create(paste0(diag_dir,'Factors/'))

Plot_factors( Report=Report, 
              ParHat=Obj$env$parList(), 
              Data=TmbData, 
              SD=Opt$SD, 
              mapdetails_list=MapDetails_List, 
              Year_Set=Year_Set, 
              category_names=levels(Data_Geostat$spp), 
              plotdir=paste0(diag_dir,'Factors/') )

############################################
## Direction of "geometric anisotropy"
## We can visualize which direction has faster or slower decorrelation 
## (termed "geometric anisotropy")
############################################
plot_anisotropy( FileName=paste0(diag_dir,"Aniso.png"), 
                 Report=Report, 
                 TmbData=TmbData )

######################################
##
######################################
