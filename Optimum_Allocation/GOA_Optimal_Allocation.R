#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################

rm(list = ls())
library(VAST); 
library(mvtnorm); library(SamplingStrata); library(sp)
library(RColorBrewer); library(raster)

VAST_model = "6c"
setwd('/Users/zackoyafuso/Google Drive/VAST_Runs/')
#setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

wd = paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
            "Optimum_Allocation/model_", VAST_model)

if(!dir.exists(wd)) dir.create(wd)

load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))

load('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/Extrapolation_depths.RData')
#load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

df = cbind(
  data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                          breaks = c(-171, -159, -154, -147, -140, -130), 
                          labels = c('Shumagin_1', 'Chirikof_2', 'Kodiak_3',
                                     'Yakutak_4', 'SE_5')),
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X = Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean)
)
names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth"),#, 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.2),
                       mut_change = c(0.01, 0.1, 0.5),
                       elitism_rate = c(0.1, 0.2, 0.5),
                       dom1 = 2:5,
                       dom2 = 2:5,
                       dom3 = 2:5,
                       dom4 = 2:5,
                       dom5 = 2.5)


ns = Save$TmbData$n_c
domains = unique(df$Domain)
ndom = length(unique(frame$domainvalue))

# rm(list  = ls()[!ls() %in% c('settings', 'frame', 'VAST_model', 'ns') ])

# for(i in 1:nrow(settings)){
#   
#   wd = paste0("C:/Users/Zack Oyafuso/Documents/",
#               "GitHub/MS_OM_GoA/Optimum_Allocation/",
#               "model_", VAST_model, "/",
#               'cv_', settings$cv[i], '_', 
#               'mut_change_', settings$mut_change[i], '_',
#               'elitism_rate_', settings$elitism_rate[i], '.RData')
i=1
cv = list()
for(spp in 1:ns) cv[[paste0('CV', spp)]] = rep(settings$cv[i], ndom)
cv[['DOM']] = levels(domains)
cv[['domainvalue']] = as.numeric(domains)
cv <- as.data.frame(cv)

set.seed(1234 + i)
solution <- optimStrata(method = "continuous",
                        errors = cv, 
                        framesamp = frame,
                        iter = 50,
                        pops = 100,
                        elitism_rate = settings$elitism_rate[i],
                        mut_chance = settings$mut_change[i],
                        nStrata = rep(5, ndom),
                        showPlot = T,
                        parallel = T)

strataStructure <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress=FALSE)

#   save(list=c('strataStructure', 'solution'), file = wd)
# }

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