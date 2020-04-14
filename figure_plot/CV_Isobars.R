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
output_dir = paste0(c('/Users/zackoyafuso/', 
                      'C:/Users/Zack Oyafuso/')[which_machine],
                    'Google Drive/MS_Optimizations/figure_plot/')

load(paste0(github_dir, 'optimization_results.RData'))

#############################
## CV Isobars
#############################
{
   png(filename = paste0(output_dir, 'CV_Isobars.png'), units = 'in',
       width = 6, height = 5, res = 500)
   par(mar = c(5,5,1,1), mfrow = c(1,1))
   plot(n ~ nstrata, data = settings, type = 'n', las = 1, 
        xlim = c(5,max(settings$nstrata) + 10), axes = F, ann = F )
   axis(side = 1, at = c(seq(5,25,5), seq(30,60,10)))
   axis(side = 2, las =1 )
   mtext(side = 1, 'Number of Strata', line = 3); mtext(side = 2, 'Total Sample Size', line = 3)
   box()
   abline(h = c(280, 550, 800), col='grey', lwd = 2, lty = 'dotted')
   
   for(icv in unique(settings$cv)[1:22 %% 2 == 1] ) {
      lines(n ~ nstrata, data = settings, subset = cv == icv)
      points(n ~ nstrata, data = settings, subset = cv == icv, 
             pch = 16, cex = 0.75)
      text(x = max(settings$nstrata[settings$cv == icv]),
           y = settings$n[settings$cv == icv][which.max(settings$nstrata[settings$cv == icv])],
           paste0(icv*100, '% CV'),
           pos = 4, cex = 0.75)
   }

   text(x = 65, y = c(300, 570, 820), c('1 Boat', '2 Boats', '3 Boats'), col = 'grey')
   #abline(lty = 'dotted', v = c(7, seq(5,25,5), seq(30,60,10)))
   dev.off()
}


#############################
## N_CV_Tradeoff
#############################
{
   png(filename = paste0(output_dir, 'N_CV_Tradeoff.png'), units = 'in',
       width = 6, height = 5, res = 500)
   par(mfrow = c(1,1), mar = c(5,5,1,1))
   plot(n ~ cv, data = settings, subset = nstrata == 5, type = 'l',
        xlab = 'Upper CV Constraint', ylab = 'Total Sample Size', las = 1)
   points(n ~ cv, data = settings, subset = (nstrata == 5),  pch = 16)
   abline(h = c(820, 550, 280), col = 'red' )
   
   opt_sol = sapply(X = c(820, 550, 280), 
                    FUN = function(x) 
                       which.min( abs(settings$n[settings$nstrata == 5]-x) ))
   
   points(n ~ cv, data = settings[(settings$nstrata == 5),][opt_sol,], col = 'red', pch = 16, cex = 2)
   legend('topright', legend = '5 strata', bty = 'n', cex = 2)
   dev.off()
}