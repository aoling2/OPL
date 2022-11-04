sim_hetero_brf <- function(i_sim, field, exp_design, coef_data) {

  result_all <- data.table()
  correlation_rho <- correlation_rho_ls[1]
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  #' ### Data preparation
  # /*~~~~~~~~~~~~~~~~~~~~~~*/
  # === merge field data with the coefs data ===#
  data <-
    coef_data[sim == i_sim, ][data.table(field), on = "cell_id"] %>%
    .[, opt_N := (- pN / pCorn + b1) / (2 * b2)] 
  
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
  
  # # === merge N treatments to field data ===#
  data_with_N <- 
    exp_design$trial_design[[1]][data, on = "plot_id"] %>% 
    N_data[., on = "n_id"] %>% 
    .[, N := N * (1 + N_error * 0.1)] %>%
    .[N < 0, N := 0] %>%
    .[, N2 := N^2]
  
  #/*----------------------------------*/
  #' ## generate yield
  #/*----------------------------------*/
  # loop over different yield monitor accuracy
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

  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  # #' ### Aggregate data by analysis unit
  # # /*~~~~~~~~~~~~~~~~~~~~~~*/
  subplot_xy <-
    data_with_y %>%
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
    .[, c('subplot_id', 'opt_N_rf_perfect', 'pi_hat')]

  for (p in 1:nrow(price_eval)) {
      pc <- price_eval[p][[1]]
      pn <- price_eval[p][[2]]
  #   # /*----------------------------------*/
  #   #' ## Economic analysis
  #   # /*----------------------------------*/

      data_return <- data_with_y %>%
        data_rf_perfect[., on = "subplot_id"] %>%
        # === true profit ===#
        .[, opt_N := (- pn / pc + b1) / (2 * b2)] %>%
        .[opt_N > Nk, opt_N := Nk] %>% 
        # === true yield and profit ===#
        .[, yield_opt := gen_yield_QP(b0, b1, b2, Nk, opt_N)] %>%
        .[, pi_opt := pc * yield_opt - pn * opt_N - fixed_c] %>%
        .[, yield_brf := gen_yield_QP(b0, b1, b2, Nk, opt_N_rf_perfect)] %>%
        .[, pi_brf := pc * yield_brf - pn * opt_N_rf_perfect - fixed_c]
      data_return[pi_opt<pi_brf,]
  #   # take field average
      new_entry <- data_return %>%
        .[, .(
          yield_accuracy = correlation_rho,
          pCorn = pc,
          pN = pn,
          opt_N_rf_perfect = mean(opt_N_rf_perfect, na.rm = T),
          pi_brf = mean(pi_brf, na.rm = T),
          opt_N = mean(opt_N, na.rm = T),
          pi_opt = mean(pi_opt, na.rm = T)
        )] %>%
        setnames(c("yield_accuracy", "pCorn", "pN", "opt_N_rf_perfect", "pi_rf", "opt_N", "pi_opt"))
      result_all <- rbind(result_all, new_entry)
    }
 return(result_all)
}