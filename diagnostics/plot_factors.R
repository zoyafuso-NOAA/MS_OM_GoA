plot_factors = function (Report, ParHat, Data, SD = NULL, Year_Set = NULL, category_names = NULL, 
    RotationMethod = "PCA", mapdetails_list = NULL, Dim_year = NULL, 
    Dim_species = NULL, plotdir = paste0(getwd(), "/"), land_color = "grey", 
    ...) 
{
    if (all(c("Options", "Options_vec") %in% names(Data))) {
        Options_vec = Data$Options_vec
        Options = Data$Options
    }
    if ("Options_list" %in% names(Data)) {
        Options_vec = Data$Options_list$Options_vec
        Options = Data$Options_list$Options
    }
    if (is.vector(Data[["FieldConfig"]]) && length(Data[["FieldConfig"]]) == 
        4) {
        Data[["FieldConfig"]] = rbind(matrix(Data[["FieldConfig"]], 
            ncol = 2, dimnames = list(c("Omega", "Epsilon"), 
                c("Component_1", "Component_2"))), Beta = c(Beta1 = -2, 
            Beta2 = -2))
    }
    else {
        if (!is.matrix(Data[["FieldConfig"]]) || !all(dim(Data[["FieldConfig"]]) == 
            c(3, 2))) {
            stop("`FieldConfig` has the wrong dimensions in `Summarize_Covariance`")
        }
    }
    if ("D_xct" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_xct)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_xct)[2]
    }
    if ("D_xcy" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_xcy)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_xcy)[2]
    }
    if ("D_gcy" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_gcy)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_gcy)[2]
    }
    Dim = function(num) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))))
    Dim_year = Dim(length(Year_Set))
    Dim_species = Dim(length(category_names))
    Psi2prime_list = Psiprime_list = Lprime_SE_list = Hinv_list = L_SE_list = Lprime_list = L_list = vector("list", 
        length = 6)
    for (i in 1:6) {
        Par_name = c("Omega1", "Epsilon1", "Beta1", "Omega2", 
            "Epsilon2", "Beta2")[i]
        if (Par_name == "Omega1") {
            Var_name = "Omegainput1_sf"
            Var2_name = "Omegainput1_gf"
        }
        if (Par_name == "Epsilon1") {
            Var_name = "Epsiloninput1_sft"
            Var2_name = "Epsiloninput1_gft"
        }
        if (Par_name == "Beta1") {
            Var_name = "beta1_ft"
            Var2_name = "missing"
        }
        if (Par_name == "Omega2") {
            Var_name = "Omegainput2_sf"
            Var2_name = "Omegainput2_gf"
        }
        if (Par_name == "Epsilon2") {
            Var_name = "Epsiloninput2_sft"
            Var2_name = "Epsiloninput2_gft"
        }
        if (Par_name == "Beta2") {
            Var_name = "beta2_ft"
            Var2_name = "missing"
        }
        if (as.vector(Data[["FieldConfig"]])[i] > 0) {
            L_list[[i]] = calc_cov(L_z = ParHat[[paste0("L_", 
                tolower(Par_name), "_z")]], n_f = as.vector(Data[["FieldConfig"]])[i], 
                n_c = Data$n_c, returntype = "loadings_matrix")
            rownames(L_list[[i]]) = category_names
            if (class(SD) == "sdreport") {
                rowindex = grep(paste0("L_", tolower(Par_name), 
                  "_z"), rownames(SD$cov.fixed))
                L_rz = mvtnorm::rmvnorm(n = 1000, mean = ParHat[[paste0("L_", 
                  tolower(Par_name), "_z")]], sigma = SD$cov.fixed[rowindex, 
                  rowindex])
                L_rcf = array(NA, dim = c(nrow(L_rz), dim(L_list[[i]])))
                for (rI in 1:nrow(L_rz)) {
                  L_rcf[rI, , ] = calc_cov(L_z = L_rz[rI, ], 
                    n_f = as.vector(Data[["FieldConfig"]])[i], 
                    n_c = Data$n_c, returntype = "loadings_matrix")
                }
                Lmean_cf = apply(L_rcf, MARGIN = 2:3, FUN = mean)
                Lsd_cf = apply(L_rcf, MARGIN = 2:3, FUN = sd)
                L_SE_list[[i]] = Lsd_cf
                rownames(L_SE_list[[i]]) = category_names
            }
            Psi_sjt = ParHat[[Var_name]]
            Psi_gjt = Report[[Var2_name]]
            if (Var_name %in% c("beta1_ft", "beta2_ft")) {
                Psi_sjt <- t(Psi_sjt)
            }
            if (is.null(Psi_sjt)) {
                stop(paste("Covariance is empty for parameter", 
                  Var_name))
            }
            logkappa = unlist(ParHat[c("logkappa1", "logkappa2")])[c(1, 
                1, 1, 2, 2, 2)[i]]
            if (Options_vec[8] == 0) {
                tau = 1/(exp(logkappa) * sqrt(4 * pi))
            }
            else if (Options_vec[8] == 1) {
                tau = 1/sqrt(1 - exp(logkappa * 2))
            }
            else stop("Check 'Options_vec[8]' for allowable entries")
            Var_rot = rotate_factors(L_pj = L_list[[i]], Psi = Psi_sjt/tau, 
                RotationMethod = RotationMethod, testcutoff = 1e-04)
            Report_tmp = list(D_xct = Var_rot$Psi_rot, Epsilon1_sct = Var_rot$Psi_rot, 
                Epsilon2_sct = Var_rot$Psi_rot)
            Lprime_list[[i]] = Var_rot$L_pj_rot
            rownames(Lprime_list[[i]]) = category_names
            Psiprime_list[[i]] = Var_rot$Psi_rot
            Hinv_list[[i]] = Var_rot$Hinv
            if (class(SD) == "sdreport") {
                rowindex = grep(paste0("L_", tolower(Par_name), 
                  "_z"), rownames(SD$cov.fixed))
                L_rz = mvtnorm::rmvnorm(n = 1000, mean = ParHat[[paste0("L_", 
                  tolower(Par_name), "_z")]], sigma = SD$cov.fixed[rowindex, 
                  rowindex])
                Lprime_rcf = array(NA, dim = c(nrow(L_rz), dim(L_list[[i]])))
                for (rI in 1:nrow(L_rz)) {
                  tmpmat = calc_cov(L_z = L_rz[rI, ], n_f = as.vector(Data[["FieldConfig"]])[i], 
                    n_c = Data$n_c, returntype = "loadings_matrix")
                  Lprime_rcf[rI, , ] = rotate_factors(L_pj = tmpmat, 
                    RotationMethod = "PCA", testcutoff = 1e-04, 
                    quiet = TRUE)$L_pj_rot
                }
                Lmean_cf = apply(Lprime_rcf, MARGIN = 2:3, FUN = mean)
                Lsd_cf = apply(Lprime_rcf, MARGIN = 2:3, FUN = sd)
                Lprime_SE_list[[i]] = Lsd_cf
                rownames(Lprime_SE_list[[i]]) = category_names
            }
            if (!is.null(Psi_gjt)) {
                Var2_rot = rotate_factors(L_pj = L_list[[i]], 
                  Psi = Psi_gjt/tau, RotationMethod = RotationMethod, 
                  testcutoff = 1e-04)
                Report_tmp = list(D_xct = Var2_rot$Psi_rot, Epsilon1_sct = Var2_rot$Psi_rot, 
                  Epsilon2_sct = Var2_rot$Psi_rot)
                Psi2prime_list[[i]] = Var2_rot$Psi_rot
            }
            Dim_factor = Dim(as.vector(Data[["FieldConfig"]])[i])
            png(file = paste0(plotdir, "Factor_loadings--", Par_name, 
                ".png"), width = Dim_factor[2] * 4, height = Dim_factor[1] * 
                4, units = "in", res = 200)
            par(mfrow = Dim_factor, mar = c(0, 2, 2, 0))
            for (cI in 1:as.vector(Data[["FieldConfig"]])[i]) FishStatsUtils::plot_loadings(L_pj = Var_rot$L_pj_rot, 
                whichfactor = cI)
            dev.off()
            if (!is.null(mapdetails_list)) {
                if (Par_name %in% c("Epsilon1", "Epsilon2")) {
                  plot_maps(plot_set = c(6, 6, NA, 7, 7, NA)[i], 
                    Report = Report_tmp, PlotDF = mapdetails_list[["PlotDF"]], 
                    MapSizeRatio = mapdetails_list[["MapSizeRatio"]], 
                    working_dir = plotdir, Year_Set = Year_Set, 
					Years2Include = Years2Include,
                    category_names = paste0("Factor_", 1:dim(Var_rot$Psi_rot)[2]), 
                    legend_x = mapdetails_list[["Legend"]]$x/100, 
                    legend_y = mapdetails_list[["Legend"]]$y/100, 
                    ...)
                }
                if (Par_name %in% c("Omega1", "Omega2")) {
                  plot_variable(Y_gt = array(Report_tmp$D_xct[, 
                    , 1], dim = dim(Report_tmp$D_xct)[1:2]), 
                    map_list = mapdetails_list, working_dir = plotdir, 
                    panel_labels = paste0("Factor_", 1:dim(Var_rot$Psi_rot)[2]), 
                    file_name = paste0("Factor_maps--", Par_name))
                }
            }
        }
        else {
            Lprime_SE_list[[i]] = L_SE_list[[i]] = L_SE_list[[i]] = Psi2prime_list[[i]] = Psiprime_list[[i]] = Lprime_list[[i]] = L_list[[i]] = "Element not estimated, and therefore empty"
        }
    }
    names(Hinv_list) = names(Psi2prime_list) = names(Psiprime_list) = names(Lprime_SE_list) = names(L_SE_list) = names(Lprime_list) = names(L_list) = c("Omega1", 
        "Epsilon1", "Beta1", "Omega2", "Epsilon2", "Beta2")
    Return = list(Loadings = L_list, Rotated_loadings = Lprime_list, 
        Rotated_factors = Psiprime_list, Rotated_projected_factors = Psi2prime_list, 
        Rotation_matrices = Hinv_list)
    if (class(SD) == "sdreport") {
        Return[["Loadings_SE"]] = L_SE_list
        Return[["Rotated_loadings_SE"]] = Lprime_SE_list
    }
    return(invisible(Return))
}
