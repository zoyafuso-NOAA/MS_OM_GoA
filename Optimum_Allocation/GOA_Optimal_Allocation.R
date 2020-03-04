
library(VAST); library(mvtnorm); library(SamplingStrata)

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')
# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')

VAST_model = "1b"
load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

Ests = Report$Index_gcyl[,,Years2Include,]

depth_01 = Extrapolation_List$Data_Extrap$depth

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
                         nStrata = nstrata3,
                         suggestions = initial_solution3)

strataStructure <- summaryStrata(solution3$framenew,
                                 solution3$aggr_strata,
                                 progress=FALSE)
strataStructure

colors = c('black', 'red', 'green', 'blue', 'orange', 'gold', 'brown')
par(mar = c(0,0,0,0))
plot(Spatial_List$loc_g, pch = 16, cex = 0.5, 
     col = colors[solution3$indices[,'X1']])
expected_CV(strata = solution3$aggr_strata)
