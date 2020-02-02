setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/Optimization_BeringSea')

library(glpkAPI)

library(tidyr); library(reshape2)

EBS = read.csv('EBS_data_trimmed.csv')

spp = length(unique(EBS$SPECIES_CODE))

spp_mean = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = mean), key = SPECIES_CODE, value = WGTCPUE)

spp_sd = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = sd), key = SPECIES_CODE, value = WGTCPUE)

#############################
## Create a function that inputs a vector of species weights
## Calculates the weighted return and total variance across spp. for a station
#############################

calc_portfolio = function(weights = rep(1/5, 5)){
  
  port_ret = matrix(nrow = nrow(spp_mean), ncol = 2,
                    dimnames = list(spp_mean$STATIONID, 
                                    c('return', 'variance')))
  
  for(id in spp_mean$STATIONID){
    temp = spread(data = subset(x = EBS, subset = STATIONID == id), 
                  key = SPECIES_CODE,
                  value = WGTCPUE)
    
    scaled_mean = colMeans(temp[,-c(1:2)]) / apply(X = spp_mean[,-1], 
                                                   MARGIN = 2, 
                                                   FUN = max)
    
    scaled_sd = apply(temp[,-c(1:2)], 
                      MARGIN = 2, 
                      FUN = sd) / apply(X = spp_sd[,-1], 
                                        MARGIN = 2, 
                                        FUN = max)
    
    temp_cor = suppressWarnings(cor(x = temp[,-c(1:2)]))
    temp_cor[is.na(temp_cor)] = 0
    
    port_var = matrix(nrow = spp, ncol = spp)
    for(i in 1:spp){
      for(j in 1:spp){
        port_var[i,j] = weights[i] * weights[j] * scaled_sd[i] * scaled_sd[j] * temp_cor[i,j]
      }
    }
    
    port_ret[id, ] = c('return' = sum(weights * scaled_mean), 
                       'var' = sum(port_var) )
  }
  
  return(port_ret)
}


do_optim = function(objvals = optim_df$return,
                    variances = optim_df$variance,
                    number_of_stations = 100,
                    var_constraint = 0.1){
  
  ##############################
  ## Defining the Model
  ###############################
  
  # Initialize Model
  model <- glpkAPI::initProbGLPK()
  
  # Set objective function as a minimization funtion
  glpkAPI::setObjDirGLPK(lp = model, 
                         lpdir = glpkAPI::GLP_MAX)
  
  # Initialize decision variables (columns)
  glpkAPI::addColsGLPK(lp = model, 
                       ncols = length(objvals))
  
  # Set the objective function, specify no bounds on decision variables
  # GLP_FR means free variable
  glpkAPI::setColsBndsObjCoefsGLPK(lp = model, 
                                   j = seq_along(objvals),
                                   lb = NULL, ub = NULL,
                                   obj_coef = objvals,
                                   type = rep(glpkAPI::GLP_FR, length(objvals)))
  
  # Specify that decision variables are binary. GLP_BV means binary variable
  glpkAPI::setColsKindGLPK(lp = model, 
                           j = seq_along(objvals),
                           kind = rep(glpkAPI::GLP_BV, length(objvals)))
  
  # Initialize the structural constraints (rows)
  # There is 1 constraint for the total variance, another for the number of stations
  glpkAPI::addRowsGLPK(lp = model, 
                       nrows = 2)
  
  mr = rep(1:2, each = length(objvals))
  mc = rep(1:length(objvals), times = 2)
  mz = c(variances, rep(1, length(objvals)))
  
  # set non-zero elements of constraint matrix
  glpkAPI::loadMatrixGLPK(lp = model, 
                          ne = length(mz),
                          ia = mr, ja = mc, ra = mz)
  
  # Set the lower and upper bounds for the right-hand side
  # of the inequalities
  lower = c(sum(variances)*var_constraint, number_of_stations)
  
  upper = c(Inf, number_of_stations)
  
  # Specify the type of structural variable. Integers refer to different types.   
  # See ?glpkAPI::glpkConstants() for the set of variable types. 
  # "2" refers to a variable with a lower bound and 
  # "3" refers a variable with an upper bound
  # "5" refers to a fixed variable
  
  bound_type = c(2, #total variance has a lower bound
                 5) #Fixed variable
  
  glpkAPI::setRowsBndsGLPK(lp = model, 
                           i = seq_along(upper),
                           lb = lower, ub = upper,
                           type = bound_type )
  
  # Presolve and automatically calculate relaxed solution
  # otherwise glpkAPI::solveSimplexGLPK(model) must be called first
  glpkAPI::setMIPParmGLPK(parm = PRESOLVE, val = GLP_ON)
  glpkAPI::setMIPParmGLPK(parm = MSG_LEV, val = GLP_MSG_ALL)
  
  # Set the maximum optimality gap of the solution
  glpkAPI::setMIPParmGLPK(parm = MIP_GAP, val = 0.001)
  
  # Stop after specified number of seconds, convert to milliseconds
  glpkAPI::setMIPParmGLPK(parm = TM_LIM, val = 1000 * 60) 
  
  #############################
  ## Solve Model
  #############################
  screen_out = glpkAPI::solveMIPGLPK(lp = model) 
  
  #############################
  ## Prepare return object
  #############################
  results <- list(output_code = screen_out,
                  status = glpkAPI::return_codeGLPK(screen_out),
                  objval = glpkAPI::mipObjValGLPK(model),
                  x = glpkAPI::mipColsValGLPK(model) )
  
  results$tot_var = sum(variances[results$x == 1])
  results$rel_var = results$tot_var / sum(variances)
  
  #Delete model
  glpkAPI::delProbGLPK(lp = model)
  
  return(results)
}

optim_df = calc_portfolio(weights = c(0.2, 0.2, 0.2, 0.2, 10) / sum(c(0.2, 0.2, 0.2, 0.2, 10)))

stations = read.csv('10110_stations_.csv')
stations = stations[stations$STATIONID != 'J-13',]
optim_df = cbind(optim_df, stations[,c('lat', 'long')])

par(mfrow = c(2,1), mar = c(0,0,0,0))
plot(lat ~ long, data = optim_df, cex = return * 20, axes = F); box()
plot(lat ~ long, data = optim_df, cex = variance * 100, axes = F); box()

res_df = res_mat = data.frame()

for(i in seq(from=50, to=300, by=25)){
  temp = 0.1; opt_res = 0
  
  while(opt_res %in% c(0, 14)){
    x = do_optim(number_of_stations = i,
                 var_constraint = temp)
    
    opt_res = x$output_code
    
    if(x$output_code %in% c(0, 14)){
      res_df = rbind(res_df, data.frame(n = sum(x$x == 1),
                                        tot_var = x$tot_var,
                                        rel_var = x$rel_var,
                                        tot_mean = x$objval) )
      
      res_mat = rbind(res_mat, as.integer(x$x))
      
      temp = x$rel_var + 0.01
      
    }
    
  }
}

res_mat = as.matrix(res_mat)

dev.off()
plot(tot_mean ~ tot_var, data = res_df, type = 'n', las = 1,
     xlim = range(res_df$tot_var), ylim = range(res_df$tot_mean),
     xlab = 'Total Variance', ylab = "Total Mean")

for(i in seq(from=50, to=300, by=25)){
  lines(tot_mean ~ tot_var, data = res_df, subset = n == i, lwd = 2)
  points(tot_mean ~ tot_var, data = res_df, subset = n == i, pch = 16)
  with(subset(res_df, n == i),
       text(max(tot_var), min(tot_mean), paste('n =', i), pos = 1)
  )
}

###########################
## Show solutions for a particular effort level
###########################
temp_df = subset(res_df, n == 100)
temp_mat = matrix( res_mat[res_df$n == 100,], ncol = nrow(optim_df) )

par(mfrow = c(3,3), mar = c(4,4,1,1))
plot(tot_mean ~ tot_var, data = temp_df, 
     type = 'b', las = 1, pch = 16, lwd = 2, cex = 2, 
     col = c('black', 'red', 'darkgreen', 'blue', 'cyan', 'magenta', 'orange', 'grey'),
     xlim = range(temp_df$tot_var), ylim = range(temp_df$tot_mean),
     xlab = 'Total Variance', ylab = "Total Mean")

par(mar = rep(0.5, 4))
for(i in 1:nrow(temp_mat)){
  plot(lat ~ long, data = stations, ann = F, axes = F, asp = 1); box()
  points(lat ~ long, data = stations[temp_mat[i,] == 1,], pch = 16, cex = 2,
         col = c('black', 'red', 'darkgreen', 'blue', 
                 'cyan', 'magenta', 'orange', 'grey')[i])
}

################################
## Species Plots
################################
dev.off()
subset(x = res_df, subset = (n == 200))
temp_mat = res_mat[res_df$n == 200,]

ispp = colnames(spp_mean)[2]

plot(x = colSums(spp_mean[temp_mat[i,] == 1,-1]) / colSums(spp_mean[,-1]),
     y = colSums(spp_var[temp_mat[i,] == 1,-1])/ colSums(spp_var[,-1]), 
     type = 'n', las = 1, 
     xlim = c(0, sum(spp_mean[,ispp])), 
     ylim = c(0,sum(spp_var[,ispp])),
     xlab = 'Mean', ylab = 'Variance')

for(i in 1:nrow(temp_mat)){
  points(x = sum(spp_mean[temp_mat[i,] == 1, ispp]),
         y = sum(spp_var[temp_mat[i,] == 1, ispp]),
         col = c('black', 'red', 'darkgreen', 'blue', 
                 'cyan', 'magenta', 'orange', 'grey')[i],
         pch = 16)
}

