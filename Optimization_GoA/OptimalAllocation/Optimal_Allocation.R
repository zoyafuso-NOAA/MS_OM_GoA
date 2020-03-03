
library(VAST); library(mvtnorm)

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/')

VAST_model = "3a"
load(paste0('VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0('VAST_output',VAST_model,'/Spatial_Settings.RData'))

source('diagnostics/plot_factors.R')
source('diagnostics/plot_residuals.R')
source('diagnostics/summarize_covariance.r')

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj


Ests = array(data = Report$Index_gcyl[attributes(Opt$SD$value)$names == 'Index_gcyl'], 
             dim =c(dim(Report$D_gcy), dim(Report$Index_cyl)[3]) )

str(TmbData$X_gtp)

depth_01 = (TmbData$X_gtp[,1,1] - min(TmbData$X_gtp[,1,1])) / diff(range(TmbData$X_gtp[,1,1]))

df = cbind(
  data.frame(Domain = "GoA",
             x = 1:250,
             depth = depth_01),
  Ests[,,1,1]
)

names(df)[-(1:3)] = paste0('spp', 1:10)

frame1 <- buildFrameDF(df = df,
                       id = "x",
                       X = c("depth"),
                       Y = paste0('spp', 1:10),
                       domainvalue = "Domain")
strata1 <- buildStrataDF(frame1, progress=F)

cv = list()
for(i in 1:10) cv[[paste0('CV', i)]] = 0.4
cv[['DOM']] = 'GoA'; cv[['domainvalue']]=1
cv <- as.data.frame(cv)
cv

checkInput(errors = checkInput(errors = cv, 
                               strata = strata1, 
                               sampframe = frame1))

allocation <- bethel(strata1,cv[1,])
sum(allocation)

set.seed(1234)
solution1 <- optimStrata(method = "atomic",
                         errors = cv, 
                         nStrata = 10,
                         framesamp = frame1,
                         iter = 50,
                         pops = 10)
expected_CV(solution1$aggr_strata)


frame2 <- buildFrameDF(df = df,
                       id = "x",
                       X = c("x"),
                       Y = paste0('spp', 1:10),
                       domainvalue = "Domain")
head(frame2)
strata2 <- buildStrataDF(frame2, progress=F)
initial_solution2 <- KmeansSolution(strata = strata2,
                                    errors = cv,
                                    maxclusters = 10)  

nstrata2 <- tapply(initial_solution2$suggestions,
                   initial_solution2$domainvalue,
                   FUN=function(x) length(unique(x)))
nstrata2
set.seed(1234)
solution2 <- optimStrata(method = "atomic",
                         errors = cv, 
                         framesamp = frame2,
                         iter = 50,
                         pops = 10,
                         nStrata = nstrata2,
                         suggestions = initial_solution2)
expected_CV(solution2$aggr_strata)

frame3 <- buildFrameDF(df = df,
                       id = "x",
                       X = "depth",
                       Y = paste0('spp', 1:10),
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
head(strataStructure)

eval3 <- evalSolution(frame = solution3$framenew, 
                      outstrata = solution3$aggr_strata, 
                      nsampl = 200,
                      progress = FALSE) 

eval3$coeff_var

eval2 <- evalSolution(frame = solution2$framenew, 
                      outstrata = solution2$aggr_strata, 
                      nsampl = 200,
                      progress = FALSE) 

eval2$coeff_var

eval1 <- evalSolution(frame = solution1$framenew, 
                      outstrata = solution1$aggr_strata, 
                      nsampl = 200,
                      progress = FALSE) 

eval3$coeff_var

colors = c('black', 'red', 'green', 'blue', 'orange', 'gold', 'brown')
plot(Spatial_List$loc_g, pch = 16, cex = 1, 
     col = colors[solution3$indices[,'X1']])
plot(Spatial_List$loc_g, pch = 16, cex = 1, 
     col = colors[solution2$indices])
plot(Spatial_List$loc_g, pch = 16, cex = 1, 
     col = colors[solution1$indices])
