extract_value = function(Sdreport, Report, Obj, variable_name, 
                         plot_value = "estimate", n_samples) {
  if (missing(Report)) {
    Report = Obj$report()
  }
  if (is.function(plot_value)) {
    if (missing(Obj)) 
      stop("Must provide `Obj` for `extract_value(.)` in `plot_maps(.)` when specifying a function for argument `plot_value`")
    Var_r = sample_variable(Sdreport = Sdreport, Obj = Obj, 
                            variable_name = variable_name, n_samples = n_samples)
    Return = apply(Var_r, MARGIN = 1:(length(dim(Var_r)) - 
                                        1), FUN = plot_value)
    if (any(dim(Return) != dim(Report[[variable_name]]))) {
      stop("Check `extract_value(.)` in `plot_maps(.)`")
    }
  }
  else if (plot_value == "estimate") {
    Return = Report[[variable_name]]
  }
  else stop("Check input `plot_value` in `plot_maps(.)`")
  return(Return)
}

setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output2a')

library(VAST); library(mvtnorm)

# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')
load(paste0('VAST_MS_GoA_Run.RData'))
load(paste0('Spatial_Settings.RData'))

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj

load('VAST_MS_GoA_Run.RData')

Opt = Save$Opt
Report = Save$Report

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

eps = Report$Epsilon2_gct


# plot_maps(plot_set=c(3), 
Sdreport=Opt$SD
category_names=levels(Data_Geostat[,'spp'])

mean_crit = apply(X = log(Report$D_gcy), MARGIN = 1:2, FUN = mean)

