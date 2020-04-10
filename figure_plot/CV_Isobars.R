#########################
## CV Isobars
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine],
                    'GitHub/MS_OM_GoA/Optimum_Allocation/', 'model_', modelno, '/')
load(paste0(github_dir, 'optimization_results.RData'))

par(mar = c(5,5,1,1))
plot(n ~ nstrata, data = settings, type = 'n', las = 1, 
     xlim = c(5,max(settings$nstrata) + 5), axes = F, ann = F )
axis(side = 1, at = c(seq(5,25,5), seq(30,60,10)))
axis(side = 2, las =1 )
mtext(side = 1, 'Number of Strata', line = 3); mtext(side = 2, 'Total Sample Size', line = 3)
box()

for(icv in unique(settings$cv)[1:22 %% 2 == 1] ) {
 lines(n ~ nstrata, data = settings, subset = cv == icv)
 points(n ~ nstrata, data = settings, subset = cv == icv, pch = 16)
 text(x = max(settings$nstrata[settings$cv == icv]),
      y = settings$n[settings$cv == icv][which.max(settings$nstrata[settings$cv == icv])],
 paste0(icv*100, '% CV'),
 pos = 4, cex = 1)
}

#abline(lty = 'dotted', v = c(7, seq(5,25,5), seq(30,60,10)))

