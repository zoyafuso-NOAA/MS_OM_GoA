
dev.off(); rm(list = ls())
library(VAST); library(mvtnorm); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

VAST_model = "1b"
load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))
load("C:/Users/Zack Oyafuso/Documents/Github/MS_OM_GoA/Extrapolation_depths.RData")

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Ests = Report$Index_gcyl[,,Years2Include,]

depth_01 = Extrapolation_depths$depth

df = cbind(
  data.frame(Domain = "GoA",
             x = 1:TmbData$n_g,
             depth = depth_01),
  apply(X = Ests, MARGIN = 1:2, FUN = mean)
)

names(df)[-(1:3)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame1 <- buildFrameDF(df = df,
                       id = "x",
                       X = c("depth"),
                       Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                       domainvalue = "Domain")
strata1 <- buildStrataDF(frame1, progress=F)

cv = list()
for(i in 1:TmbData$n_c) cv[[paste0('CV', i)]] = 0.1
cv[['DOM']] = 'GoA'; cv[['domainvalue']]=1
cv <- as.data.frame(cv)
cv

checkInput(errors = checkInput(errors = cv, 
                               strata = strata1, 
                               sampframe = frame1))

allocation <- bethel(strata1,cv[1,])
sum(allocation)

frame3 <- buildFrameDF(df = df,
                       id = "x",
                       X = "depth",
                       Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                       domainvalue = "Domain")

head(frame3)
set.seed(1234)
init_sol3 <- KmeansSolution2(frame=frame3,
                             errors=cv,
                             maxclusters = 10)  

nstrata3 <- tapply(init_sol3$suggestions,
                   init_sol3$domainvalue,
                   FUN=function(x) length(unique(x)))
nstrata3

initial_solution3 <- prepareSuggestion(init_sol3,frame3,nstrata3)
set.seed(1234)
solution3 <- optimStrata(method = "continuous",
                         errors = cv, 
                         framesamp = frame3,
                         iter = 50,
                         pops = 10,
				 minnumstr=10,
                         nStrata = nstrata3,
                         suggestions = initial_solution3,
parallel = TRUE, cores = 6)

strataStructure <- summaryStrata(solution3$framenew,
                                 solution3$aggr_strata,
                                 progress=FALSE)
strataStructure

colors = c('black', 'red', 'green', 'blue', 'orange', 'gold', 'brown')
plot(Spatial_List$loc_g, pch = 16, cex = 0.5, 
     col = colors[solution3$indices[,'X1']])
expected_CV(strata = solution3$aggr_strata)


## Spatial Method
index_means =apply(X = Ests, MARGIN = 1:2, FUN = mean)
index_vars = apply(X = Ests, MARGIN = 1:2, FUN = var)
colnames(index_means) = paste0(gsub(x = Save$Spp, pattern = ' ', replacement = '_'), '_mean')
colnames(index_vars) = paste0(gsub(x = Save$Spp, pattern = ' ', replacement = '_'), '_var')

df = cbind(
  data.frame(Domain = "GoA",
		dom1 = 1,
		 depth = Extrapolation_depths$depth,
             x = 1:TmbData$n_g,
		lon = Extrapolation_depths$E_km,
		lat = Extrapolation_depths$N_km),
  index_means, index_vars
  
)

df <- as.data.frame(df)
head(df)

frame <- buildFrameSpatial(df=df,
                      id="x",
                      X="depth",
                      Y=paste0(gsub(x = Save$Spp, pattern = ' ', replacement = '_'), '_mean'),
                      variance=paste0(gsub(x = Save$Spp, pattern = ' ', replacement = '_'), '_var'),
                      lon="lon",
                      lat="lat",
                      domainvalue = "dom1")

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.1,1),
                         CV2=rep(0.1,1),
CV3=rep(0.1,1),
CV4=rep(0.1,1),
                         domainvalue=c(1:1) ))

init_sol <- KmeansSolution2(frame=frame,
                             errors=cv,
                             maxclusters = 10)  

nstrata <- tapply(init_sol$suggestions,
                   init_sol$domainvalue,
                   FUN=function(x) length(unique(x)))
nstrata
initial_solution <- prepareSuggestion(init_sol,frame,nstrata)

set.seed(1234)
solution <- optimStrata (
  method = "spatial",
  errors=cv, 
  framesamp=frame,
  iter = 15,
  pops = 10,
  nStrata = 10,
  fitting = rep(1, 4),
  range = rep(1, 4),
  kappa=1,
  writeFiles = FALSE,
  showPlot = TRUE)

plot(Extrapolation_depths[,c("Lon", 'Lat')], pch = '.', cex = 1,
col = heat.colors(4)[cut(frame$Y1, breaks = quantile(frame$Y1))] )

strataStructure <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress=FALSE)
strataStructure

colors = c('black', 'red', 'green', 'blue', 'orange', 'gold', 'brown', 'yellow')
plot(Spatial_List$loc_g, pch = 16, cex = 0.5, 
     col = colors[solution$indices[,'X1']])
expected_CV(strata = solution$aggr_strata)

save(list=ls()[!ls() %in% c('Extrapolation_List', 'gulf_of_alaska_grid', 'Ests', 'df',
					'Obj', 'Opt', 'Save', 'Report', 'Data_Geostat', 
					"Extrapolation_depths", 'Spatial_List', 'TmbData') ],
file = paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/test.RData"))