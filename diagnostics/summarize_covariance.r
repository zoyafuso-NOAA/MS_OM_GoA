# Report 
# Data = TmbData 
# ParHat = Obj$env$parList()
# SD = Opt$SD
# category_order = 1:Data$n_c
# category_names = 1:Data$n_c
# plotdir = diag_dir
# figname = "Cov"
# plotTF = NULL
# plot_cor = TRUE
# mgp = c(2, 0.5, 0)
# tck = -0.02
# oma = c(0, 5, 2, 0)

summarize_covariance = function (Report, Data, ParHat, SD = NULL,
                                 category_order = 1:Data$n_c,
                                 category_names = 1:Data$n_c,
                                 plotdir = paste0(getwd(), "/"),
                                 figname = "Cov",
                                 plotTF = NULL,
                                 plot_cor = TRUE,
                                 mgp = c(2, 0.5, 0),
                                 tck = -0.02,
                                 oma = c(0, 5, 2, 0), ...)
{
  
  if (is.vector(Data[["FieldConfig"]]) && length(Data[["FieldConfig"]]) == 
      4) {
    Data[["FieldConfig"]] = rbind(matrix(Data[["FieldConfig"]], 
                                         ncol = 2, dimnames = list(c("Omega", "Epsilon"), 
                                                                   c("Component_1", "Component_2"))), Beta = c(Beta1 = -2, 
                                                                                                               Beta2 = -2))
  } else {
    if (!is.matrix(Data[["FieldConfig"]]) || !all(dim(Data[["FieldConfig"]]) == 
                                                  c(3, 2))) {
      stop("`FieldConfig` has the wrong dimensions in `Summarize_Covariance`")
    }
  }
  if (is.null(plotTF)) {
    plotTF = as.vector(Data[["FieldConfig"]] > 0)
  } else {
    plotTF = as.vector(plotTF)
  }
  
  Return = list()
  for (i in which(Data[["FieldConfig"]] >= 0)) {
    Par_name = c("omega1", "epsilon1", "beta1", "omega2", 
                 "epsilon2", "beta2")[i]
    L_name = paste0("L_", Par_name, "_z")
    if (!is.null(SD)) {
      sd_summary = summary(SD)
      Slot_name = paste0("lowercov_uppercor_", Par_name)
      if (Slot_name %in% rownames(sd_summary)) {
        Cor = Cov = Mat = ThorsonUtilities::Extract_SE(SD = SD, 
                                                       parname = Slot_name, columns = 1:2, Dim = c(Data$n_c, 
                                                                                                   Data$n_c))
        dimnames(Cor) = dimnames(Cov) = list(category_names, 
                                             category_names, c("Estimate", "Std.Error"))
        Cor[, , 1][lower.tri(Cor[, , 1])] = t(Mat[, , 
                                                  1])[lower.tri(Mat[, , 1])]
        diag(Cor[, , 1]) = 1
        Cor[, , 2][lower.tri(Cor[, , 2])] = t(Mat[, , 
                                                  2])[lower.tri(Mat[, , 2])]
        diag(Cor[, , 2]) = NA
        Cov[, , 1][upper.tri(Cov[, , 1])] = t(Mat[, , 
                                                  1])[upper.tri(Mat[, , 1])]
        Cov[, , 2][upper.tri(Cov[, , 2])] = t(Mat[, , 
                                                  2])[upper.tri(Mat[, , 2])]
      }
      else {
        Cov = Cor = NULL
      }
    }
    else {
      Cov = Cor = NULL
    }
    if (is.null(Cov) | is.null(Cor)) {
      Cov = Cor = array(NA, dim = c(Data$n_c, Data$n_c, 
                                    2), dimnames = list(category_names, category_names, 
                                                        c("Estimate", "Std.Error")))
      Cov[, , "Estimate"] = FishStatsUtils:::calc_cov(L_z = ParHat[[L_name]], 
                                                      n_f = as.vector(Data[["FieldConfig"]])[i], n_c = Data$n_c)
      Cor[, , "Estimate"] = cov2cor(Cov[, , "Estimate"])
    }
    List = list(Cor, Cov)
    names(List) = paste0(c("Cor_", "Cov_"), Par_name)
    Return = c(Return, List)
  }
  
  if (!is.null(figname)) {
    Dim = c(3, 2)
    if (sum(ifelse(plotTF > 0, 1, 0)) == 1) 
      Dim = c(1, 1)
    if (all(ifelse(plotTF > 0, 1, 0) == c(1, 1, 0, 0, 0, 
                                          0)) | all(ifelse(plotTF > 0, 1, 0) == c(0, 0, 1, 
                                                                                  1, 0, 0))) 
      Dim = c(1, 2)
    if (all(ifelse(plotTF > 0, 1, 0) == c(1, 0, 1, 0, 0, 
                                          0)) | all(ifelse(plotTF > 0, 1, 0) == c(0, 1, 0, 
                                                                                  1, 0, 0))) 
      Dim = c(2, 1)
    if (plot_cor == TRUE) {
      convert = function(Cov) ifelse(is.na(cov2cor(Cov)), 
                                     0, cov2cor(Cov))
    }
    else {
      convert = function(Cov) ifelse(is.na(Cov), 0, Cov)
    }
    ThorsonUtilities::save_fig(file = paste0(plotdir, figname, 
                                             "--Analytic.png"), width = Dim[2] * 4 + 1, height = Dim[1] * 
                                 4)#, ...)
    par(mfrow = Dim, mar = c(0, 1, 1, 0), mgp = mgp, tck = tck, 
        oma = oma)
    for (i in 1:6) {
      if (i %in% which(plotTF > 0)) {
        Cov_cc = FishStatsUtils:::calc_cov(L_z = ParHat[c("L_omega1_z", 
                                                          "L_epsilon1_z", "L_beta1_z", "L_omega2_z", 
                                                          "L_epsilon2_z", "L_beta2_z")][[i]], n_f = as.vector(Data[["FieldConfig"]])[i], 
                                           n_c = Data$n_c)
        plot_cov(Cov = convert(Cov_cc)[category_order, 
                                       category_order], names = list(category_names[category_order], 
                                                                     NA)[[ifelse(i == 1 | i == 3 | Dim[2] == 1, 
                                                                                 1, 2)]], names2 = list(1:nrow(Cov_cc), NA)[[ifelse(i == 
                                                                                                                                      1 | i == 2, 1, 2)]], digits = 1, font = 2)
      }
    }
    dev.off()
  }
  return(invisible(Return))
  
}