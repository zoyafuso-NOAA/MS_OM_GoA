##############################
## Bethel Algorithm on GoA current strata
##############################
rm(list = ls())

############################
## Import Libraries
############################
library(rgdal); library(raster); library(rgeos); library(tidyr)
library(SamplingStrata)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI' = 3)[1]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine], 'GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/VAST_Runs/',
                    'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/',
                    'C:\\Users\\zack.oyafuso\\Desktop\\VAST_Runs\\'),
                  'VAST_output', modelno)[which_machine]

############################
## Load Data
############################
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
load(paste0(github_dir, '/Optimum_Allocation/model_', modelno, 
            '/optimization_data_model_', modelno, '.RData'))

#Survey data
survey_data = read.csv(paste0(github_dir, '/data/data/',
                              'GOA_multspp_with_strata.csv'))
survey_data = subset(survey_data, SPECIES_NAME != 'Anoplopoma fimbria')
survey_data$SPECIES_NAME = droplevels(survey_data$SPECIES_NAME)

## Constants
sci_names = levels(Data_Geostat$spp)
N = nrow(frame)

#Calculate weights of each stratum
stratapop = table(gulf_of_alaska_grid$GOA_STRATUM)

#Calculate samples allocated across strata across years
samples_by_str = with(subset(survey_data, 
                             SPECIES_NAME == 'Sebastolobus alascanus'),
                      table(YEAR, STRATUM))
strata = colnames(samples_by_str)


######################
## Results Objects
######################

CVs = seq(0.20, 0.60, by = 0.01)
sample_allocation = array(data = 0, dim = c(length(CVs), length(strata)),
                          dimnames = list(CVs, strata))

spp_cvs = array(data = 0, dim = c(length(CVs), ns),
                dimnames = list(CVs, sci_names) )

#######
sample_mean = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM,
                                      data = survey_data, FUN=mean, drop = F),
                     key = SPECIES_NAME, value = CPUE)[,-1]
names(sample_mean) = paste0('M', 1:ns)

sample_var = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM,
                                     data = survey_data, FUN = var, drop = F),
                    key = SPECIES_NAME, value = CPUE)[,-1]
names(sample_var) = paste0('S', 1:ns)

sample_mean[is.na(sample_var)] = 0
sample_var[is.na(sample_var)] = 0

df = cbind(data.frame(stratum = strata,
                      N = as.vector(stratapop)),#,
           #X1 = factor(1:length(temp_strata))),
           sample_mean, sqrt(sample_var),
           data.frame(cens = 0,
                      cost = 1,
                      DOM1 = 'tot'))

for(icv in CVs){
 stmt = paste0('cbind(', paste0("CV", 1:ns, '=', icv, collapse = ', '), ')' )
 CV = eval(parse(text = stmt))
 errors = cbind(data.frame(DOM = 'DOM1'), CV, domainvalue = 1)
 
 n = bethel(stratif = df, errors = errors, printa=TRUE, epsilon = 1e-11, maxiter = 200)
 sample_allocation[paste(icv),] = as.numeric(n)
 spp_cvs[paste(icv),] = as.numeric( attributes(n)$outcv[,'ACTUAL CV'] )
}

total_sample_size = apply(sample_allocation, MARGIN = 1, sum)

{png(filename = paste0(github_dir, '/Optimum_Allocation_AI/GOA_test_ST.png'),
     width = 8, height = 5, units = 'in', res = 500)
 par(mfrow = c(1,1), mar = c(5,5,1,1))
 plot(1, pch = 16, las = 1, xlim = c(0.20,0.60), ylim = c(0, 1600), 
      type = 'n', xlab = 'Upper CV Constraint', ylab = 'Total Sample Size')
 lines(CVs, total_sample_size)
 points(CVs, total_sample_size, pch = 16)
 
 abline(h = c(280, 550, 820), lty = 'dotted')
 text(x = .59, y = c(350, 620, 880), c('1 Boat', paste(2:3, 'Boats')))
 dev.off()}


{png(filename = paste0(github_dir, 
                       '/Optimum_Allocation_AI/GOA_test_sppCV.png'),
     width = 8, height = 7, units = 'in', res = 500)
par(mfrow = c(1,1), mar = c(4,4,1,1))
G1 = c(3,9,13,14)
matplot(CVs, spp_cvs[,G1], type = 'l', ,ylim = c(0,0.6), col = 'black',
        las = 1, xlim = c(0.2, 0.8), lty = 1, lwd = 2,
        axes = F, xlab = 'Upper CV constraint', ylab = 'Species-specific CV')
segments(x0 = 0.2, x1 = 0.6, y0 = 0, lwd = 3, col = 'grey')
points(y = rep(0,3), 
       x = CVs[sapply(X = c(820, 550, 280), 
                  FUN = function(x) 
                   which.min(abs(total_sample_size - x)))], 
       pch = 16, col = 'grey' )
segments(x0 = CVs[sapply(X = c(820, 550, 280), 
                         FUN = function(x) 
                          which.min(abs(total_sample_size - x)))],
         x1 = CVs[sapply(X = c(820, 550, 280), 
                         FUN = function(x) 
                          which.min(abs(total_sample_size - x)))],
         y0 = 0, y1 = 0.6, col = 'grey', lty = 'dotted')
text(y = rep(0,3), 
       x = CVs[sapply(X = c(820, 550, 280), 
                      FUN = function(x) 
                       which.min(abs(total_sample_size - x)))], 
       c('1 Boat', '2 Boats', '3 Boats'), pos = 1, col = 'grey' )

axis(side = 1, at = seq(0.2, 0.6, 0.1))
axis(side = 2, las = 1)
box()

segments(x0 = 0.61, y0 = min(spp_cvs[nrow(spp_cvs),G1]), 
         y1 = max(spp_cvs[nrow(spp_cvs),G1]), lwd = 3)
text(x = 0.61, y = mean(spp_cvs[nrow(spp_cvs),G1]), 
     paste0(paste0(c('Choke Species:', sci_names[G1[1:(length(G1)-1)]]), '\n', collapse = ''), sci_names[G1[length(G1)]], sep = ''),
     pos = 4, cex = 0.75)


G2 = c(11:12)
matlines(CVs, spp_cvs[,G2], col = 'red', lty= 1, lwd = 2)
#matpoints(CVs, spp_cvs[,G2], pch = 16, col = 'black')
segments(x0 = 0.61, y0 = min(spp_cvs[nrow(spp_cvs),G2]), 
         y1 = max(spp_cvs[nrow(spp_cvs),G2]), lwd = 3)
text(x = 0.61, y = mean(spp_cvs[nrow(spp_cvs),G2]), 
     paste0(paste0(c('Initial Choke Species,\nThen Asymptotes:', sci_names[G2[1:(length(G2)-1)]]), '\n', collapse = ''), sci_names[G2[length(G2)]], sep = ''),
     pos = 4, cex = 0.75)

G3 = c(2,6,8)#,10,15) #Linear Trend
matlines(CVs, spp_cvs[,G3], col = 'blue', lty= 1, lwd = 2)
#matpoints(CVs, spp_cvs[,G3], pch = 16, col = 'black')
segments(x0 = 0.61, y0 = min(spp_cvs[nrow(spp_cvs),G3]), 
         y1 = max(spp_cvs[nrow(spp_cvs),G3]), lwd = 3)
text(x = 0.61, y = mean(spp_cvs[nrow(spp_cvs),G3]), 
     paste0(paste0(c('Linear Increase:', sci_names[G3[1:(length(G3)-1)]]), '\n', collapse = ''), sci_names[G3[length(G3)]], sep = ''),
     pos = 4, cex = 0.75)

G4 = c(1,4,5,7,10,15) #small change
matlines(CVs, spp_cvs[,G4], col = 'green', lty= 1, lwd = 2)
#matpoints(CVs, spp_cvs[,G4], pch = 16, col = 'black')
segments(x0 = 0.61, y0 = min(spp_cvs[nrow(spp_cvs),G4]), 
         y1 = max(spp_cvs[nrow(spp_cvs),G4]), lwd = 3)
text(x = 0.61, 
     y = quantile(spp_cvs[nrow(spp_cvs),G4], probs = 0.15), 
     paste0(paste0(c('Gradual Increase:', sci_names[G4[1:(length(G4)-1)]]), '\n', collapse = ''), sci_names[G4[length(G4)]], sep = ''),
     pos = 4, cex = 0.75)


dev.off()}
