plot_maps_density = function (plot_set = 3, Report, PlotDF, Sdreport = NULL, TmbData = NULL, 
    projargs = "+proj=longlat", Panel = "Category", 
    Year_Set = NULL, Years2Include = NULL, category_names = NULL, 
    quiet = FALSE, working_dir = paste0(getwd(), "/"), 
    MapSizeRatio, n_cells, ...) 
{
    if ("D_xt" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:ncol(Report$D_xt)
        if (is.null(Years2Include)) 
            Years2Include = 1:ncol(Report$D_xt)
        category_names = "singlespecies"
        Ncategories = length(category_names)
        Nyears = dim(Report$D_xt)[2]
    }
    if ("D_xct" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_xct)[3]
        if (is.null(Years2Include)) 
            Years2Include = 1:dim(Report$D_xct)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_xct)[2]
        Ncategories = dim(Report$D_xct)[2]
        Nyears = dim(Report$D_xct)[3]
    }
    if ("D_xcy" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_xcy)[3]
        if (is.null(Years2Include)) 
            Years2Include = 1:dim(Report$D_xcy)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_xcy)[2]
        Ncategories = dim(Report$D_xcy)[2]
        Nyears = dim(Report$D_xcy)[3]
    }
    if ("D_gcy" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$D_gcy)[3]
        if (is.null(Years2Include)) 
            Years2Include = 1:dim(Report$D_gcy)[3]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$D_gcy)[2]
        Ncategories = dim(Report$D_gcy)[2]
        Nyears = dim(Report$D_gcy)[3]
    }
    if ("dhat_ktp" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$dhat_ktp)[2]
        if (is.null(Years2Include)) 
            Years2Include = 1:dim(Report$dhat_ktp)[2]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$dhat_ktp)[3]
        Ncategories = dim(Report$dhat_ktp)[3]
        Nyears = dim(Report$dhat_ktp)[2]
    }
    if ("dpred_ktp" %in% names(Report)) {
        if (is.null(Year_Set)) 
            Year_Set = 1:dim(Report$dpred_ktp)[2]
        if (is.null(Years2Include)) 
            Years2Include = 1:dim(Report$dpred_ktp)[2]
        if (is.null(category_names)) 
            category_names = 1:dim(Report$dpred_ktp)[3]
        Ncategories = dim(Report$dpred_ktp)[3]
        Nyears = dim(Report$dpred_ktp)[2]
    }
    if (missing(MapSizeRatio)) {
        MapSizeRatio = c(3, 3)
    }
    if (Nyears != length(Year_Set)) {
        stop("Problem with `Year_Set`")
    }
    if (Ncategories != length(category_names)) {
        stop("Problem with `category_names`")
    }
    Return = NULL
    for (plot_num in plot_set) {
        Array_xct = NULL
        plot_code <- c("encounter_prob", "pos_catch", 
            "density", "", "", "epsilon_1", 
            "epsilon_2", "linear_predictor_1", "linear_predictor_2", 
            "density_CV", "covariates", "total_density", 
            "covariate_effects_1", "covariate_effects_2", 
            "omega_1", "omega_2")[plot_num]
        if (plot_num == 1) {
            if (quiet == FALSE) 
                message(" # Plotting presence/absense maps")
            if ("D_xt" %in% names(Report)) 
                Array_xct = Report$R1_xt
            if ("D_xct" %in% names(Report)) 
                Array_xct = Report$R1_xct
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$R1_xcy
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$R1_gcy
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
            message("`plot_num=1` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells")
        }
        if (plot_num == 2) {
            if (quiet == FALSE) 
                message(" # Plotting positive catch rate maps")
            if ("D_xt" %in% names(Report)) 
                Array_xct = log(Report$R2_xt)
            if ("D_xct" %in% names(Report)) 
                Array_xct = log(Report$R2_xct)
            if ("D_xcy" %in% names(Report)) 
                Array_xct = log(Report$R2_xcy)
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$R2_gcy
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
            message("`plot_num=2` doesn't work well when using ObsModel[2]==1, because average area-swept doesn't generally match area of extrapolation-grid cells")
        }
        if (plot_num == 3) {
            if (quiet == FALSE) 
                message(" # Plotting density maps")
            if ("D_xt" %in% names(Report)) 
                Array_xct = log(Report$D_xt)
            if ("D_xct" %in% names(Report)) 
                Array_xct = log(Report$D_xct)
            if ("D_xcy" %in% names(Report)) 
                Array_xct = log(Report$D_xcy)
            if ("D_gcy" %in% names(Report)) 
                Array_xct = log(Report$D_gcy)
            if ("dhat_ktp" %in% names(Report)) 
                Array_xct = aperm(Report$dhat_ktp[, , cI], c(1, 
                  3, 2))
            if ("dpred_ktp" %in% names(Report)) 
                Array_xct = aperm(Report$dpred_ktp[, , cI], c(1, 
                  3, 2))
        }
        if (plot_num == 4) {
            stop("`plot_num=4` is deprecated")
        }
        if (plot_num == 5) {
            stop("`plot_num=5` is deprecated")
        }
        if (plot_num == 6) {
            if (quiet == FALSE) 
                message(" # Plotting spatio-temporal effects (Epsilon) in 1st linear predictor")
            if ("D_xt" %in% names(Report)) 
                Array_xct = Report$Epsilon1_st
            if ("D_xct" %in% names(Report)) 
                Array_xct = Report$Epsilon1_sct
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$Epsilon1_sct
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$Epsilon1_gct
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
        }
        if (plot_num == 7) {
            if (quiet == FALSE) 
                message(" # Plotting spatio-temporal effects (Epsilon) in 2nd linear predictor")
            if ("D_xt" %in% names(Report)) 
                Array_xct = Report$Epsilon2_st
            if ("D_xct" %in% names(Report)) 
                Array_xct = Report$Epsilon2_sct
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$Epsilon2_sct
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$Epsilon2_gct
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
        }
        if (plot_num == 8) {
            if (quiet == FALSE) 
                message(" # Plotting 1st predictor after action of link function")
            if ("D_xt" %in% names(Report)) 
                Array_xct = Report$P1_xt
            if ("D_xct" %in% names(Report)) 
                Array_xct = Report$P1_xct
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$P1_xcy
            if ("D_gcy" %in% names(Report)) 
                stop("`plot_maps` not implemented for requested plot_num")
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
        }
        if (plot_num == 9) {
            if (quiet == FALSE) 
                message(" # Plotting 2nd predictor after action of link function")
            if ("D_xt" %in% names(Report)) 
                Array_xct = Report$P2_xt
            if ("D_xct" %in% names(Report)) 
                Array_xct = Report$P2_xct
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$P2_xcy
            if ("D_gcy" %in% names(Report)) 
                stop("`plot_maps` not implemented for requested plot_num")
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("Not implemented for SpatialVAM")
        }
        if (plot_num == 10) {
            if (quiet == FALSE) 
                message(" # Plotting density maps")
            if (is.null(Sdreport)) 
                stop("Must supply 'Sdreport' if 'plot_num=10'")
            if ("D_xt" %in% names(Report)) {
                if (!("log(Index_xtl)" %in% rownames(TMB::summary.sdreport(Sdreport)))) 
                  stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'SpatialDeltaGLMM'")
                Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == 
                  "log(Index_xtl)"), ], dim = c(dim(Report$D_xt), 
                  ncol(Report$Index_tl), 2), dimnames = list(NULL, 
                  NULL, NULL, c("Estimate", "Std. Error")))[, 
                  , 1, "Std. Error"]
            }
            if ("D_xct" %in% names(Report)) {
                if (!("log(Index_xctl)" %in% rownames(TMB::summary.sdreport(Sdreport)))) 
                  stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
                Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == 
                  "log(Index_xctl)"), ], dim = c(dim(Report$D_xct), 
                  dim(Report$Index_ctl)[3], 2), dimnames = list(NULL, 
                  NULL, NULL, NULL, c("Estimate", "Std. Error")))[, 
                  , , 1, "Std. Error"]
            }
            if ("D_xcy" %in% names(Report)) {
                if (!("log(Index_xcyl)" %in% rownames(TMB::summary.sdreport(Sdreport)))) 
                  stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
                Array_xct = array(TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport)) == 
                  "log(Index_xcyl)"), ], dim = c(dim(Report$D_xcy), 
                  dim(Report$Index_cyl)[3], 2), dimnames = list(NULL, 
                  NULL, NULL, NULL, c("Estimate", "Std. Error")))[, 
                  , , 1, "Std. Error"]
            }
            if (any(c("dhat_ktp", "dpred_ktp") %in% 
                names(Report))) 
                stop("'plot_num=10' not implemented for 'SpatialVAM'")
            Array_xct = sqrt(exp(Array_xct^2) - 1)
            if ("D_gcy" %in% names(Report)) 
                stop("`plot_maps` not implemented for requested plot_num")
        }
        if (plot_num == 11) {
            if (quiet == FALSE) 
                message(" # Plotting covariates")
            if (is.null(TmbData)) 
                stop("Must provide `TmbData` to plot covariates")
            if ("X_xtp" %in% names(TmbData)) 
                Array_xct = aperm(TmbData$X_xtp, perm = c(1, 
                  3, 2))
            if ("X_gtp" %in% names(TmbData)) 
                Array_xct = aperm(TmbData$X_gtp, perm = c(1, 
                  3, 2))
            category_names = 1:dim(Array_xct)[2]
        }
        if (plot_num == 12) {
            if (quiet == FALSE) 
                message(" # Plotting total density")
            if ("D_xt" %in% names(Report)) 
                Array_xct = log(Report$D_xt)
            if ("D_xct" %in% names(Report)) 
                Array_xct = log(apply(Report$D_xct, FUN = sum, 
                  MARGIN = c(1, 3)))
            if ("D_xcy" %in% names(Report)) 
                Array_xct = log(apply(Report$D_xcy, FUN = sum, 
                  MARGIN = c(1, 3)))
            if ("D_gcy" %in% names(Report)) 
                Array_xct = log(apply(Report$D_gcy, FUN = sum, 
                  MARGIN = c(1, 3)))
            logsum = function(vec) {
                max(vec) + log(sum(exp(vec - max(vec))))
            }
            if ("dhat_ktp" %in% names(Report)) 
                Array_xct = apply(aperm(Report$dhat_ktp, c(1, 
                  3, 2)), FUN = logsum, MARGIN = c(1, 3))
            if ("dpred_ktp" %in% names(Report)) 
                Array_xct = apply(aperm(Report$dpred_ktp, c(1, 
                  3, 2)), FUN = logsum, MARGIN = c(1, 3))
        }
        if (plot_num == 13) {
            if (quiet == FALSE) 
                message(" # Plotting covariate effects for 1st linear predictor")
            if ("D_xt" %in% names(Report)) 
                stop()
            if ("D_xct" %in% names(Report)) 
                stop()
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$eta1_xct
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$eta1_gct
            if ("dhat_ktp" %in% names(Report)) 
                stop()
            if ("dpred_ktp" %in% names(Report)) 
                stop()
        }
        if (plot_num == 14) {
            if (quiet == FALSE) 
                message(" # Plotting covariate effects for 2nd linear predictor")
            if ("D_xt" %in% names(Report)) 
                stop()
            if ("D_xct" %in% names(Report)) 
                stop()
            if ("D_xcy" %in% names(Report)) 
                Array_xct = Report$eta2_xct
            if ("D_gcy" %in% names(Report)) 
                Array_xct = Report$eta2_gct
            if ("dhat_ktp" %in% names(Report)) 
                stop()
            if ("dpred_ktp" %in% names(Report)) 
                stop()
        }
        if (is.null(Array_xct)) 
            stop("Problem with `plot_num` in `plot_maps(.)")
        if (tolower(Panel) == "category") {
            if (length(dim(Array_xct)) == 2) 
                Nplot = 1
            if (length(dim(Array_xct)) == 3) 
                Nplot = dim(Array_xct)[2]
            for (cI in 1:Nplot) {
                if (length(dim(Array_xct)) == 2) 
                  Return = Mat_xt = Array_xct
                if (length(dim(Array_xct)) == 3) 
                  Return = Mat_xt = array(as.vector(Array_xct[, 
                    cI, ]), dim = dim(Array_xct)[c(1, 3)])
                file_name = paste0(plot_code, ifelse(Nplot > 
                  1, paste0("--", category_names[cI]), 
                  ""))
                plot_args = plot_variable_density(Y_gt = Mat_xt[, Years2Include, 
                  drop = FALSE], map_list = list(PlotDF = PlotDF, 
                  MapSizeRatio = MapSizeRatio), projargs = projargs, 
                  working_dir = working_dir, panel_labels = Year_Set[Years2Include], 
                  file_name = file_name, n_cells = n_cells, ...)
            }
        }
        if (tolower(Panel) == "year") {
            Nplot = length(Years2Include)
            for (tI in 1:Nplot) {
                if (length(dim(Array_xct)) == 2) 
                  Mat_xc = Array_xct[, Years2Include[tI], drop = TRUE]
                if (length(dim(Array_xct)) == 3) 
                  Mat_xc = Array_xct[, , Years2Include[tI], drop = TRUE]
                Return = Mat_xc = array(as.vector(Mat_xc), dim = c(dim(Array_xct)[1], 
                  Ncategories))
                file_name = paste0(plot_code, ifelse(Nplot > 
                  1, paste0("--", Year_Set[Years2Include][tI]), 
                  ""))
                plot_args = plot_variable_density(Y_gt = Mat_xc, map_list = list(PlotDF = PlotDF, 
                  MapSizeRatio = MapSizeRatio), projargs = projargs, 
                  working_dir = working_dir, panel_labels = category_names, 
                  file_name = file_name, n_cells = n_cells, ...)
            }
        }
    }
    if (is.null(Return) & quiet == FALSE) 
        message(" # No plots selected in `plot_set`")
    return(invisible(Return))
}