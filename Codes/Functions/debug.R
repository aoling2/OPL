i_sim <- 1
i <- 7
{
  temp_design <- trial_designs[i, ]
  field_col <- temp_design$field.col
  
  field_raw <- readRDS(paste0("Data/field_col_", field_col, ".rds"))
  
  coef_data <- 
    readRDS(
      here(
        "Data", 
        paste0(
          'field_col_400', 
          field_col,
          '_coef.rds'
        )
      )
    )
  correlation_rho <- temp_design$correlation_rho
  # /*----------------------------------*/
  #' ## Assign plot, block id
  # /*----------------------------------*/
  # === generate field id
  field <-
    gen_field_ids(
      field_config = exp_designs[i, ],
      field_raw
    )
  exp_design <- temp_design
  # /*----------------------------------*/
  #' ## Simulations
  # /*----------------------------------*/
  # === simulations ===#
}
sim <- function(i_sim, field, field_raw, field_col, exp_design, coef_data, correlation_rho) {
  
  result_all <- data.table()
  
  # if (field_col==576){
  #   
  #   # === spatial feature ===# 
  #   gstat_model <- "Sph"
  #   sp_range <- 600
  #   
  #   xy <-
  #     st_centroid(field_raw) %>%
  #     st_coordinates() %>%
  #     data.table() %>%
  #     .[, cell_id := field_raw$cell_id]
  #   
  #   num_sim_p <- 1
  #   
  #   coef_data <- para_gen(xy = xy, 
  #                         sp_range = sp_range, 
  #                         gstat_model = gstat_model, 
  #                         num_sim_p = num_sim_p)
  #   
  #   data <-
  #     coef_data[data.table(field), on = "cell_id"] %>%
  #     .[, opt_N := (- pN / pCorn + b1) / (2 * b2)]
  #   #=== homo ===#
  #   set.seed(847259)
  #   r <- sample(1:nrow(field), 1)
  #   data_homo <- data.table(field) %>%
  #     .[, b0:=coef_data[sim==i_sim,][r,]$b0] %>%
  #     .[, b1:=coef_data[sim==i_sim,][r,]$b1] %>%
  #     .[, b2:=coef_data[sim==i_sim,][r,]$b2] %>%
  #     .[, Nk:=coef_data[sim==i_sim,][r,]$Nk] %>%
  #     .[, plateau:=coef_data[sim==i_sim,][r,]$plateau]%>%
  #     .[, m_error:=coef_data[sim==i_sim,]$m_error] %>%
  #     .[, N_error:=coef_data[sim==i_sim,]$N_error] 
  #   
  # }else{
    # /*~~~~~~~~~~~~~~~~~~~~~~*/
    #' ### Data preparation
    # /*~~~~~~~~~~~~~~~~~~~~~~*/
    # === merge field data with the coefs data ===#
    data <-
      coef_data[sim == i_sim, ][data.table(field), on = "cell_id"] %>%
      .[, opt_N := (- pN / pCorn + b1) / (2 * b2)] 
    # /*~~~~~~~~~~~~~~~~~~~~~~*/
    #' Homo data 
    # /*~~~~~~~~~~~~~~~~~~~~~~*/
    #=== pick a random number r  ===#
    set.seed(8472259)
    r <- sample(1:nrow(field), 1)
    #=== merge field data with the coefs data ===#
    data_homo <- data.table(field) %>%
      .[, b0:=coef_data[sim==i_sim,][r,]$b0] %>%
      .[, b1:=coef_data[sim==i_sim,][r,]$b1] %>%
      .[, b2:=coef_data[sim==i_sim,][r,]$b2] %>%
      .[, Nk:=coef_data[sim==i_sim,][r,]$Nk] %>%
      .[, plateau:=coef_data[sim==i_sim,][r,]$plateau]%>%
      .[, m_error:=coef_data[sim==i_sim,]$m_error] %>%
      .[, N_error:=coef_data[sim==i_sim,]$N_error] 
  # }
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### Assign N
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  N_levels <- quantile(
    data[, opt_N] %>% pmin(., 300),
    prob = seq(0.05, 0.95, length = exp_design$num_treatment)
  )
  
  N_data <-
    N_levels%>%
    data.table(N = .) %>%
    .[, n_id := 1:nrow(.)]
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' Homo N treatment
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  N_levels_homo <- if (exp_design$num_treatment==4) {
    N_homo_4} else if (exp_design$num_treatment==6) {
      N_homo_6
    } else {
      N_homo_8
    }
  N_data_homo <- 
    N_levels_homo%>% 
    data.table(N = .) %>% 
    .[, n_id := 1:nrow(.)]
  
  # # === merge N treatments to field data ===#
  data_with_N <-
    exp_design$trial_design[[1]][data, on = "plot_id"] %>%
    N_data[., on = "n_id"] %>%
    .[, N := N * (1 + N_error * 0.1)] %>%
    .[N < 0, N := 0] %>%
    .[, N2 := N^2]
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' Homo merge N
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  data_with_N_homo <- 
    exp_design$trial_design[[1]][data_homo, on = "plot_id"] %>% 
    N_data_homo[., on = "n_id"] %>% 
    .[, N := N * (1 + N_error * 0.1)] %>%
    .[N < 0, N := 0] %>%
    .[, N2 := N^2]
  
  #/*----------------------------------*/
  #' ## generate yield
  #/*----------------------------------*/
  conv_factor <- (1 - correlation_rho^2) / correlation_rho^2
  
  data_with_y <-
    data_with_N[, det_yield := gen_yield_QP(b0, b1, b2, Nk, N)] %>%
    # === add error component ===#
    .[, yerror := det_yield * m_error] %>%
    .[, yield_sp := det_yield * (1 + m_error)] %>%
    .[, machineyerror := rnorm(nrow(data), 0, sqrt(conv_factor * var(yield_sp)))] %>%
    .[, yield := yield_sp + machineyerror] %>%
    # === keep the relevant vars ===#
    .[, .(
      cell_id, plot_id, subplot_id, N,
      N2, yield, b0, b1, b2,
      Nk, plateau, opt_N, type,
      X, Y, geometry
    )]
  
  data_with_y_homo <- 
    data_with_N_homo[, det_yield := gen_yield_QP(b0, b1, b2, Nk, N)] %>%
    # === add error component ===#
    .[, yerror := det_yield * m_error] %>%
    .[, yield_sp := det_yield * (1 + m_error)] %>%
    # .[, machineyerror := rnorm(nrow(data), 0, sqrt(conv_factor * var(yield_sp)))] %>%
    .[, machineyerror := rnorm(nrow(data), 0, 50*sqrt(conv_factor * var(yield_sp)))] %>%
    .[, yield := yield_sp + machineyerror] %>%
    # === keep the relevant vars ===#
    .[, .(
      cell_id, plot_id, subplot_id, N,
      N2, yield, b0, b1, b2,
      Nk, plateau, type,
      X, Y, geometry
    )]
  
  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  # #' ### Clean Data by 4 sd
  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  # # data_clean <- clean_yield(data_with_y, "yield")
  
  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  # #' ### Aggregate data by analysis unit
  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  subplot_xy <-
    data_with_y_homo %>%
    .[, .(X = mean(X), Y = mean(Y)), by = subplot_id]
  
  # # # === by analysis unit ===#
  reg_data <-
    data_with_y[
      type == "Analysis",
      .(
        yield = mean(yield),
        N = mean(N),
        N2 = mean(N2),
        b0 = mean(b0),
        b1 = mean(b1),
        b2 = mean(b2),
        Nk = mean(Nk)
      ),
      by = subplot_id
    ] %>%
    left_join(subplot_xy, by = "subplot_id")
  
  reg_data_homo <-
    data_with_y_homo[
      type == "Analysis",
      .(
        yield = mean(yield),
        N = mean(N),
        N2 = mean(N2),
        b0 = mean(b0),
        b1 = mean(b1),
        b2 = mean(b2),
        Nk = mean(Nk)
      ),
      by = subplot_id
    ] %>%
    left_join(subplot_xy, by = "subplot_id")
  
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### GWR
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  gwr_beta <- GWR(reg_data)
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### BRF
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  
  #=== use perfect info ===#
  X <- reg_data %>% data.table() %>% .[, c("N","b0","b1","b2","Nk"), with = FALSE] %>% data.frame()
  Y <- reg_data %>% data.table() %>% .[, yield]
  BRF_temp <- boosted_regression_forest(
    X = X,
    Y = Y,
    num.trees = 2000,
    min.node.size = 10)
  N_seq <- seq(min(N_levels), max(N_levels), by = 2)
  data_test <- reg_data %>% data.table()
  data_rf_perfect <-
    data_test[, c('subplot_id',"b0","b1","b2","Nk", 'yield'), with = FALSE] %>%
    .[rep(1:nrow(.), each = length(N_seq)), ] %>%
    .[, N := rep(N_seq, nrow(.)/length(N_seq))] %>%
    .[, yield_hat := predict(BRF_temp, newdata = .[, c("N","b0","b1","b2","Nk"), with = FALSE])] %>%
    .[, pi_hat := pCorn * yield_hat - pN * N - fixed_c] %>%
    .[, .SD[which.max(pi_hat)], by = subplot_id] %>%
    .[, opt_N_rf_perfect := N] %>%
    .[, c('subplot_id', 'opt_N_rf_perfect')]
  
  #/*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### GAM
  #/*~~~~~~~~~~~~~~~~~~~~~~*/
  
  gam_res <- gam(yield~s(N, k = ifelse (exp_design$num_treatment==4, exp_design$num_treatment-1, 5) ) + 
                   s(X, k = 6) + 
                   s(Y, k = 6) + 
                   ti(X, Y, k = 6), 
                 data=reg_data_homo)
  
  ## loop over different price combinations
  for (p in 1:nrow(price_eval)) {
    pc <- price_eval[p][[1]]
    pn <- price_eval[p][[2]]
    #   # /*----------------------------------*/
    #   #' ## Economic analysis
    #   # /*----------------------------------*/
    # calculate optimal N for homo field
    opt_N_gam <- data.table(N=seq(min(N_levels_homo)-50,max(N_levels_homo),length=1000)) %>%
      .[, X := reg_data_homo[1, ] %>%  pull(X)] %>%
      .[, Y := reg_data_homo[1, ] %>%  pull(Y)] %>%
      .[, yhat := predict(gam_res,newdata=.)] %>%
      .[, profit := pc*yhat-pn*N] %>%
      .[profit == max(profit), N] %>%
      max(min(N_levels_homo),.) %>%
      min(max(N_levels_homo),.)
    
    data_return <- data_with_y %>%
      gwr_beta[., on = "subplot_id"] %>%
      data_rf_perfect[., on = "subplot_id"] %>%
      # === opt_N_gwr ===#
      .[, opt_N_gwr := (pn / pc - b1_hat) / (2 * b2_hat)] %>%
      # === limit the range of opt_N_gwr ===#
      .[opt_N_gwr < min(N_levels), opt_N_gwr := min(N_levels)] %>%
      .[opt_N_gwr > max(N_levels), opt_N_gwr := max(N_levels)] %>%
      # === yield and profit ===#
      .[, yield_gwr := gen_yield_QP(b0, b1, b2, Nk, opt_N_gwr)] %>%
      .[, pi_gwr := pc * yield_gwr - pn * opt_N_gwr - fixed_c] %>%
      # === brf yield and profit ===#
      .[, yield_brf := gen_yield_QP(b0, b1, b2, Nk, opt_N_rf_perfect)] %>%
      .[, pi_brf := pc * yield_brf - pn * opt_N_rf_perfect - fixed_c] %>%
      # === true yield and profit ===#
      .[, opt_N := (- pn / pc + b1) / (2 * b2)] %>%
      .[opt_N > Nk, opt_N := Nk] %>%
      # === yield and profit ===#
      .[, yield_opt := gen_yield_QP(b0, b1, b2, Nk, opt_N)] %>%
      .[, pi_opt := pc * yield_opt - pn * opt_N - fixed_c]
    
    data_return_homo <- data_with_y_homo %>%
      # === yield and profit ===#
      .[, yield_gam := gen_yield_QP(b0, b1, b2, Nk, opt_N_gam)] %>%
      .[, pi_opt_gam := pc * yield_gam - pn * opt_N_gam - fixed_c] %>%
      # === true yield and profit ===#
      .[, opt_N_homo := (- pn / pc + b1) / (2 * b2)] %>%
      # === yield and profit ===#
      .[, yield_opt_homo := gen_yield_QP(b0, b1, b2, Nk, opt_N_homo)] %>%
      .[, pi_opt_homo := pc * yield_opt_homo - pn * opt_N_homo - fixed_c]
    
    ## take field average
    new_entry <- data_return_homo[data_return, on = "cell_id"] %>%
      .[, .(
        i_sim = i_sim,
        yield_accuracy = correlation_rho,
        pCorn = pc,
        pN = pn,
        opt_N_gwr = mean(opt_N_gwr, na.rm = T),
        pi_gwr = mean(pi_gwr, na.rm = T),
        opt_N_rf_perfect = mean(opt_N_rf_perfect, na.rm = T),
        pi_brf = mean(pi_brf, na.rm = T),
        opt_N = mean(opt_N, na.rm = T),
        pi_opt = mean(pi_opt, na.rm = T),
        opt_N_gam=opt_N_gam,
        pi_opt_gam=mean(pi_opt_gam, na.rm=T),
        opt_N_homo=mean(opt_N_homo, na.rm=T),
        pi_opt_homo=mean(pi_opt_homo, na.rm=T)
      )] %>%
      setnames(c("i_sim", "yield_accuracy", "pCorn", "pN", "opt_N_gwr", "pi_gwr","opt_N_rf_perfect", "pi_brf", "opt_N", "pi_opt",
                 'opt_N_gam', 'pi_opt_gam', 'opt_N_homo', 'pi_opt_homo'))
    result_all <- rbind(result_all, new_entry)
  }
  return(result_all)
}
