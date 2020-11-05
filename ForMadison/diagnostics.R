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

##################################################
#### Set up directories     
##################################################
rm(list = ls())
which_machine <- c("Zack_PC" = 1, "Zack_GI_PC" = 2)[2]

github_dir <- c("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/")[which_machine]
VAST_dir <- c("C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Single_Species/",
              "G:/Oyafuso/VAST_Runs_EFH/Single_Species/")[which_machine]

depth_in_model <- T
for (ispp in c(1:2)) {
  which_spp <- c( "Sebastes polyspinis", "Sebastes variabilis", 
                  "Sebastes brevispinis", "Microstomus pacificus", 
                  "Lepidopsetta polyxystra", "Lepidopsetta bilineata",
                  "Hippoglossus stenolepis", "Hippoglossoides elassodon",
                  "Glyptocephalus zachirus", "Gadus macrocephalus",
                  "Gadus chalcogrammus", "Sebastes alutus",
                  "Atheresthes stomias", "Sebastolobus alascanus")[ispp]
  

  
  result_dir <- paste0(VAST_dir, which_spp, 
                       ifelse(depth_in_model,  "_depth", ""), "/")
  
  fun_dir <- paste0(github_dir, "diagnostics/")
  
  if(!dir.exists(paste0(result_dir, "diagnostics/")) ) 
    dir.create(paste0(result_dir, "diagnostics/"))
  
  ##################################################
  ####  Load VAST fit and import customized plot functions 
  ##################################################
  load(paste0(result_dir, "/fit.RData"))
  
  for (ifile in c("plot_residuals.R")) {
    source( paste0(fun_dir, ifile) )
    rm(ifile)
  }
  
  ##################################################
  #### Extract objects from fitted object   
  #### Set up constants
  ##################################################
  Data_Geostat <- fit$data_frame
  names(Data_Geostat)[c(1:2, 5:7)] <- c("Lat", "Lon", "Catch_KG","Year", "spp")
  
  Spatial_List <- fit$spatial_list
  Extrapolation_List <- fit$extrapolation_list
  Opt <- fit$parameter_estimates
  Obj <- fit$tmb_list$Obj
  Report <- fit$Report
  TmbData <- fit$data_list
  
  Year_Set <- seq(min(Data_Geostat[, "Year"]),max(Data_Geostat[, "Year"]))
  Years2Include <- which( Year_Set %in% sort(unique(Data_Geostat[, "Year"])))
  
  Index_Ests <- fit$Report$Index_cyl[,Years2Include,1]
  Index_SDs <- matrix(data = fit$parameter_estimates$SD$sd[attributes(fit$parameter_estimates$SD$value)$names == "Index_cyl"], nrow = 1 )[,Years2Include]
  
  
  MapDetails_List <- make_map_info( "Region" = "Gulf_of_Alaska", 
                                    "spatial_list" = Spatial_List, 
                                    "Extrapolation_List" = Extrapolation_List)
  
  xrange <- range(Extrapolation_List$Data_Extrap[, "E_km"])
  yrange <- range(Extrapolation_List$Data_Extrap[, "N_km"])
  xrange_diff <- diff(xrange)
  yrange_diff <- diff(yrange)
  
  AK <- rgdal::readOGR(paste0(dirname(fun_dir), 
                              "/data/shapefiles/AKland.shp")) 
  AK <- sp::spTransform(AK, 
                        CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                     " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                     " +ellps=GRS80 +datum=NAD83 +units=km",
                                     " +no_defs"))
  
  CA <- rgdal::readOGR(paste0(dirname(fun_dir), 
                              "/data/shapefiles/canada_dcw.shp")) 
  CA <- sp::spTransform(CA, 
                        CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                     " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                     " +ellps=GRS80 +datum=NAD83 +units=km",
                                     " +no_defs"))
  
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
  
  plot_data(Extrapolation_List = Extrapolation_List,
            Spatial_List = Spatial_List,
            Data_Geostat = Data_Geostat,
            Year_Set = Years2Include,
            PlotDir = paste0(result_dir, "diagnostics/"))
  
  ############################################
  ## Convergence
  ## Here I print the diagnostics generated during parameter estimation, and I 
  ## confirm that (1) no parameter is hitting an upper or lower bound and (2) the
  ## final gradient for each fixed-effect is close to zero. For explanation of 
  ## parameters, please see `?make_data`.
  ############################################
  # pander::pandoc.table(Opt$diagnostics[,c("Param", "Lower", "MLE",
  # "Upper", "final_gradient")] ) 
  
  # max(abs(Opt$diagnostics$final_gradient))
  
  ############################################
  ## Diagnostics for encounter-probability component
  ## Next, we check whether observed encounter frequencies for either low or high 
  ## probability samples are within the 95% predictive interval for predicted 
  ## encounter probability
  ############################################
  Enc_prob <- plot_encounter_diagnostic(Report = Report, 
                                        Data_Geostat = Data_Geostat, 
                                        DirName = paste0(result_dir, 
                                                         "diagnostics/"))
  
  ############################################
  ## Diagnostics for positive-catch-rate component
  ## We can visualize fit to residuals of catch-rates given encounters using a
  ## Q-Q plot.  A good Q-Q plot will have residuals along the one-to-one line.  
  ############################################
  
  Q <- plot_quantile_diagnostic( 
    TmbData = TmbData,
    Report = Report,
    FileName_PP = "Posterior_Predictive",
    FileName_Phist = "Posterior_Predictive-Histogram",
    FileName_QQ = "Q-Q_plot",
    FileName_Qhist = "Q-Q_hist",
    DateFile = paste0(result_dir, "diagnostics") )
  
  ############################################
  ## We then plot Pearson residuals.  If there are visible patterns (areas with 
  ## consistently positive or negative residuals accross or within years) then 
  ## this is an indication of the model "overshrinking" results towards the 
  ## intercept, and model results should then be treated with caution.  
  ############################################
  
  PResid <- plot_residuals(Lat_i = Data_Geostat[, "Lat"], 
                          Lon_i = Data_Geostat[, "Lon"], 
                          TmbData = TmbData, 
                          Report = Report, 
                          Q = Q, 
                          spatial_list = Spatial_List,
                          extrapolation_list = Extrapolation_List,
                          working_dir = paste0(result_dir, "diagnostics/"), 
                          Year_Set = Year_Set, 
                          Years2Include = Years2Include, 
                          mar = c(0, 0, 2, 0), 
                          oma = c(3.5, 3.5, 0, 0), 
                          cex = 1.8)
  
  save(list = c("Q", "PResid"),
       file = paste0(result_dir, "diagnostics/", "diagnostics.RData"))
  
  ############################################
  ## Plot Density across years for each Species
  ############################################
  {
    png(paste0(result_dir, "diagnostics/density.png"), 
        units = "in", 
        height = 5, 
        width = 7, 
        res = 500)
    
    #Plot layout
    par(mar = c(0, 0, 0, 0), 
        oma = rep(0.5, 4), 
        mfrow = c(4, 3))
    
    for(iyear in Years2Include){
      
      #Extract density values for a species in a year, 
      vals  = Report$D_gcy[,1,iyear]
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
    dev.off()
    
  }
  
  ############################################
  ## Plot Omega of the 1st and second components
  ############################################
  {
    png(filename = paste0(result_dir, "diagnostics/omega.png"),
        units = "in", height = 6, width = 12, res = 500)
    
    #Plot layout
    par(mar = c(0,0,0,0), oma = c(0,0,2,0))
    layout(mat = matrix(c(1:2), byrow = T, nrow = 1),
           widths = c(1,0.2))
    offset = 0
    
    #Empty plot
    plot(1, type = "n", axes = F,
         xlim = xrange, ylim = yrange + c(-0.5*yrange_diff, 0))
    for(omegatype in 1:2){ #Two types for the 0/1 and pos components
      
      #Extract spatial component
      scaled_var = list(Report$Omega1_gc[,1],
                        Report$Omega2_gc[,1])[[omegatype]]
      
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
    rm(offset, omegatype, scaled_var)
  }
  
  ############################################
  ## Plot Epsilon for the first, middle and last year
  ## for each of the 1st and second components
  ############################################
  {
    png(paste0(result_dir, "diagnostics/epsilon.png"),
        units = "in", height = 5, width = 6, res = 500)
    
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
        scaled_var = list(Report$Epsilon1_gct[,1,iyear],
                          Report$Epsilon2_gct[,1,iyear])[[epstype]]
        
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
    
    #Title
    mtext(side = 3, outer = T, 
          paste0("Spatiotemporal Effect in Occurrence (Top)", 
                 " and Positive Response (Bottom)"), line = 1)
    
    dev.off()
  }
  
  ############################################
  ## Index of abundance
  ## Compare VAST estimates with Design-Based (Stratified RS) estimates
  ############################################
  Index = plot_biomass_index( DirName=paste0(result_dir, "diagnostics/"), 
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
  
  {
    png(paste0(result_dir, "diagnostics/indices_comparison_DBE.png"), 
        width = 12, height = 6, units = "in", res = 500)
    par(mfrow = c(1,1), mar = c(3,3,2,1), oma = c(0,2.5,0,0))
    
    GOA_DBE = readRDS(file = paste0(github_dir,
                                    "data/GOA_biomass_indices_wnames.rds") )
    
    #Design based Estimator and SD Interval
    temp_DBE = subset(GOA_DBE, SPECIES_NAME == which_spp &
                        YEAR %in% Year_Set[Years2Include])
    
    temp_DBE = temp_DBE[order(temp_DBE$YEAR),]
    
    upper_DBE = (temp_DBE$TOTAL_BIOMASS + sqrt(temp_DBE$BIOMASS_VAR))/1e6
    lower_DBE = (temp_DBE$TOTAL_BIOMASS - sqrt(temp_DBE$BIOMASS_VAR))/1e6
    
    #VAST SD Intervals
    upper = (Index_Ests + Index_SDs)/1e6
    lower = (Index_Ests - Index_SDs)/1e6
    
    #Empty plot
    plot(x = Year_Set[Years2Include], Index_Ests/1e6, type = "n", 
         ylab = "Index", xlab = "Year", las = 1, 
         ylim = c(0, max(c(upper, upper_DBE) )) )
    
    #Plot VAST interals
    polygon(x = c(Year_Set[Years2Include],
                  rev(Year_Set[Years2Include])),
            y = c(lower, rev(upper)), col = "grey", lty = "dotted")
    
    lines(x = Year_Set[Years2Include], Index_Ests/1e6)
    points(x = Year_Set[Years2Include], Index_Ests/1e6, pch= 16)
    
    #Plot design-based intervals
    points(Year_Set[Years2Include],
           temp_DBE$TOTAL_BIOMASS/1e6, col = "red", pch = 16)
    lines(Year_Set[Years2Include],
          temp_DBE$TOTAL_BIOMASS/1e6, col = "red")
    
    segments(x0 = Year_Set[Years2Include],
             x1 = Year_Set[Years2Include],
             y0 = lower_DBE,
             y1 = upper_DBE, 
             col = "red")
    
    #Plot Legend
    # plot(1, type = "n", axes = F, ann = F)
    legend("topleft", legend = c("DBE", "VAST"), 
           col = c("red", "black"), pch = 16, lty=1, cex = 1)
    mtext(side = 2, outer = T, text = "Abundance Index (million metric tons)",
          line = 1)
    dev.off()
  }
  
  ############################################
  ## Direction of "geometric anisotropy"
  ## We can visualize which direction has faster or slower decorrelation 
  ## (termed "geometric anisotropy")
  ############################################
  plot_anisotropy( FileName=paste0(result_dir, "diagnostics/Aniso.png"), 
                   Report=Report, 
                   TmbData=TmbData )
}

