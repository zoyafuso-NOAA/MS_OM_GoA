plot_residuals = function (Lat_i, Lon_i, TmbData, Report, Q, projargs = "+proj=longlat", 
    working_dir = paste0(getwd(), "/"), spatial_list, extrapolation_list, 
    Year_Set = NULL, Years2Include = NULL, zrange, ...) 
{
    if (!("t_iz" %in% names(TmbData))) {
        TmbData$t_iz = matrix(TmbData$t_i, ncol = 1)
    }
    if (!("t_yz" %in% names(TmbData))) {
        TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol = 1)
    }
    exp_rate_xy = obs_rate_xy = total_num_xy = exp_num_xy = obs_num_xy = matrix(NA, 
        nrow = spatial_list$n_x, ncol = nrow(TmbData$t_yz))
    for (yI in 1:nrow(TmbData$t_yz)) {
        which_i_in_y = (TmbData$t_iz == outer(rep(1, TmbData$n_i), 
            TmbData$t_yz[yI, ]))
        which_i_in_y = which(apply(which_i_in_y, MARGIN = 1, 
            FUN = all))
        if (length(which_i_in_y) > 0) {
            exp_rate_xy[, yI] = tapply(Report$R1_i[which_i_in_y], 
                INDEX = factor(spatial_list$knot_i[which_i_in_y], 
                  levels = 1:spatial_list$n_x), FUN = mean)
            obs_rate_xy[, yI] = tapply(TmbData$b_i[which_i_in_y] > 
                0, INDEX = factor(spatial_list$knot_i[which_i_in_y], 
                levels = 1:spatial_list$n_x), FUN = mean)
            total_num_xy[, yI] = tapply(TmbData$b_i[which_i_in_y], 
                INDEX = factor(spatial_list$knot_i[which_i_in_y], 
                  levels = 1:spatial_list$n_x), FUN = length)
        }
        else {
            total_num_xy[, yI] = 0
        }
        exp_num_xy = exp_rate_xy * total_num_xy
        obs_num_xy = obs_rate_xy * total_num_xy
    }
    Q1_xy = (obs_num_xy - exp_num_xy)/sqrt(exp_num_xy * (total_num_xy - 
        exp_num_xy)/total_num_xy)
    which_pos = which(TmbData$b_i > 0)
    bvar_ipos = bpred_ipos = NULL
    if (all(c("var_y", "pred_y") %in% names(Q))) {
        bvar_ipos = Q[["var_y"]]
        bpred_ipos = Q[["pred_y"]]
    }
    if (all(c("var_y", "pred_y") %in% names(Q[[1]]))) {
        bvar_ipos = bpred_ipos = rep(NA, length = length(which_pos))
        for (i_e in 1:length(Q)) {
            which_pos_and_e = which(TmbData$e_i[which_pos] == 
                (i_e - 1))
            bvar_ipos[which_pos_and_e] = Q[[i_e]][["var_y"]]
            bpred_ipos[which_pos_and_e] = Q[[i_e]][["pred_y"]]
        }
    }
    if (is.null(bvar_ipos) & is.null(bpred_ipos)) {
        stop("Something is wrong with `Q` input")
    }
    sum_obs_xy = sum_exp_xy = var_exp_xy = matrix(NA, nrow = spatial_list$n_x, 
        ncol = nrow(TmbData$t_yz))
    for (yI in 1:nrow(TmbData$t_yz)) {
        which_i_in_y = (TmbData$t_iz == outer(rep(1, TmbData$n_i), 
            TmbData$t_yz[yI, ]))
        which_i_in_y = which(apply(which_i_in_y, MARGIN = 1, 
            FUN = all))
        which_i_in_y_and_pos = intersect(which_i_in_y, which_pos)
        which_ipos_in_y = (TmbData$t_iz[which_pos, ] == outer(rep(1, 
            length(which_pos)), TmbData$t_yz[yI, ]))
        which_ipos_in_y = which(apply(which_ipos_in_y, MARGIN = 1, 
            FUN = all))
        if (length(which_i_in_y_and_pos) > 0) {
            sum_obs_xy[, yI] = tapply(TmbData$b_i[which_i_in_y_and_pos], 
                INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], 
                  levels = 1:spatial_list$n_x), FUN = sum)
            sum_exp_xy[, yI] = tapply(bpred_ipos[which_ipos_in_y], 
                INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], 
                  levels = 1:spatial_list$n_x), FUN = sum)
            var_exp_xy[, yI] = tapply(bvar_ipos[which_ipos_in_y], 
                INDEX = factor(spatial_list$knot_i[which_i_in_y_and_pos], 
                  levels = 1:spatial_list$n_x), FUN = sum)
        }
    }
    Q2_xy = (sum_obs_xy - sum_exp_xy)/sqrt(var_exp_xy)
    if (!is.null(working_dir)) {
        for (zI in 1:2) {
            Q_xy = list(Q1_xy, Q2_xy)[[zI]]
            if (!missing(zrange)) {
                Q_xy = ifelse(Q_xy < zrange[1], zrange[1], Q_xy)
                Q_xy = ifelse(Q_xy > zrange[2], zrange[2], Q_xy)
                zlim = zrange
            }
            else {
                zlim = c(-1, 1) * ceiling(max(abs(Q_xy), na.rm = TRUE))
            }
            Col = colorRampPalette(colors = c("blue", "white", 
                "red"))
            textmargin = "Pearson residual"
            plot_code = c("pearson_residuals_1", "pearson_residuals_2")[zI]
            x2i = spatial_list$NN_Extrap$nn.idx[, 1]
            Include = extrapolation_list[["Area_km2_x"]] > 0 & 
                extrapolation_list[["a_el"]][, 1] > 0
            DF = cbind(extrapolation_list$Data_Extrap[, c("Lon", 
                "Lat")], x2i = x2i, Include = Include)
            if (is.null(Year_Set)) 
                Year_Set = 1:ncol(Q_xy)
            if (is.null(Years2Include)) 
                Years2Include = 1:ncol(Q_xy)
            plot_args = plot_variable(
				#Y_gt = ifelse(is.na(Q_xy), mean(zlim), Q_xy),
				Y_gt = Q_xy[,Years2Include],			
				map_list = list(PlotDF = DF), 
				legend_y = c(0.05,0.45), legend_x = c(0.7, 0.75),
                projargs = projargs, working_dir = working_dir, 
                #panel_labels = Year_Set, 
				panel_labels = Year_Set[Years2Include],
				file_name = plot_code, 
                zlim = zlim, col = Col, ...)
        }
    }
    Return = list(Q1_xy = Q1_xy, Q2_xy = Q2_xy)
    return(invisible(Return))
}


