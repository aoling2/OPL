FitFlextableToPage <- function(ft, pgwidth = 6.5){
  
  ft_out <- ft 
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

### =======================
### Generate parameters
### =======================
#/*----------------------------------*/
#' ## Generate coefs spatial field based on 'gstat' model
#/*----------------------------------*/
gen_coefs <- function(xy, mean, psill, range, coef_name, nsim){
  
  g_N <- gstat(formula=z~1, locations=~X+Y, dummy=T,
               beta=mean, model=vgm(psill=psill, range=range, nugget=0.1, 
                                    model=gstat_model), nmax=50)
  
  b_sim <- predict(g_N, newdata=xy, nsim=nsim) %>%
    data.table() %>%
    melt(id.vars=c('X','Y')) %>%
    setnames(c('variable','value'),c('sim',coef_name)) %>%
    .[,sim:=as.numeric(gsub('sim','',sim))] %>%
    xy[.,on=c('X','Y')] %>%
    .[,c("cell_id","sim","X","Y",coef_name),with=FALSE]
  
  return(b_sim)
}
#/*----------------------------------*/
#' ## Generate error
#/*----------------------------------*/
#=== use the 'Sph' model; 'Gau' is too smooth for errors
gen_errors <- function(xy, mean, psill, range, coef_name, nsim){
	
	g_N <- gstat(formula=z~1, locations=~X+Y, dummy=T,
				 beta=mean, model=vgm(psill=psill, range=range, nugget=0, 
				 					  model=gstat_model), nmax=50)
	
	b_sim <- predict(g_N, newdata=xy, nsim=nsim) %>%
		data.table() %>%
		melt(id.vars=c('X','Y')) %>%
		setnames(c('variable','value'),c('sim',coef_name)) %>%
		.[,sim:=as.numeric(gsub('sim','',sim))] %>%
		xy[.,on=c('X','Y')] %>%
		.[,c("cell_id","sim","X","Y",coef_name),with=FALSE]
	
	return(b_sim)
}
para_gen <- function(xy, sp_range, gstat_model, num_sim_p){
  # y = -b2 * x2 + b1 * x - b0
  # === N_star ===#
  N_star <- gen_coefs(xy, mean = 200, psill = 1000, range = sp_range, coef_name = "N_star", nsim = num_sim_p) %>%
    # >>> normalize <<<
    .[,sd_b:=sd(N_star),by=sim] %>%
    .[,mean_b:=mean(N_star),by=sim] %>%
    .[,p:=pnorm(N_star,mean=mean_b, sd=sd_b)] %>%
    .[,N_star:=100+p*150] %>%
    .[, c("cell_id", "sim", "N_star")]
  # === ymax ===#
  ymax <- gen_coefs(xy, mean = 15000, psill = 3000000, range = sp_range, coef_name = "ymax", nsim = num_sim_p) %>%
    # >>> normalize <<<
    .[,sd_b:=sd(ymax),by=sim] %>%
    .[,mean_b:=mean(ymax),by=sim] %>%
    .[,p:=pnorm(ymax,mean=mean_b, sd=sd_b)] %>%
    .[,ymax:=8000+p*8000] %>%
    .[, c("cell_id", "sim", "ymax")]
  # === b0 ===#
  b0 <- gen_coefs(xy, mean = 6000, psill = 200000, range = sp_range, coef_name = "b0", nsim = num_sim_p) %>%
    .[ymax, on = c("sim", "cell_id")] %>%
    # >>> normalize <<< (b0 needs to be < ymax)
    .[,sd_b:=sd(b0),by=sim] %>%
    .[,mean_b:=mean(b0),by=sim] %>%
    .[,p:=pnorm(b0,mean=mean_b, sd=sd_b)] %>%
    .[,b0:=3000+p*4000] %>%
    .[, c("cell_id", "sim", "b0")]
  # === b1, b2 ===#
  coef_data <- N_star %>%
    .[ymax, on = c("sim", "cell_id")] %>%
    .[b0, on = c("sim", "cell_id")] %>%
    # === derive b1, b2 from b0, ymax, and N_star
    .[, b2 := (ymax - b0) / N_star^2] %>%
    .[, b1 := 2 * N_star * b2] 

  # === plateau/Nk ===# (simply use N_star as Nk)
  coef_data <- coef_data %>%
    .[, Nk := N_star] %>%
    .[, plateau := ymax]%>%
    .[, c("cell_id", "sim", "b0", "b1", "b2", "Nk", "plateau")]
  # ====m errors====#
  m_error <- gen_errors(xy, mean = 0, psill = 0.015, range = sp_range, coef_name = "m_error", nsim = num_sim_p) %>%
    # >>> normalize <<<
    .[, min_b := min(m_error), by = sim] %>%
    .[, max_b := max(m_error), by = sim] %>%
    .[, p := punif(m_error, min_b, max_b)] %>%
    .[, m_error := -0.3 + p * 0.6] %>%
    .[, c("cell_id", "sim", "m_error")]
  # ====cell-level N errors====#
  N_error <- gen_errors(xy, mean = 0, psill = 0.2, range = 50, coef_name = "N_error", nsim = num_sim_p) %>%
    # >>> normalize <<<
    .[, min_b := min(N_error), by = sim] %>%
    .[, max_b := max(N_error), by = sim] %>%
    .[, p := punif(N_error, min_b, max_b)] %>%
    .[, N_error := -1 + p * 2] %>%
    .[, c("cell_id", "sim", "N_error")]
  # === combine all ===#
  par_data <- m_error[coef_data, on = c("cell_id", "sim")]%>%.[N_error, on = c("cell_id", "sim")]
  return(par_data)
}
#/*=================================================*/
#' # Quadratic-Plateau response
#/*=================================================*/
gen_yield_QP <- function(b_0,b_1,b_2,Nk,N){
  yield <- (N<Nk)*(b_0+b_1*N-b_2*N^2) + (N>=Nk)*(b_0+b_1*Nk-b_2*Nk^2)
  return(yield)
}
#/*=================================================*/
#' # Creating data frame for graph
#/*=================================================*/
expand_grid_df <- function(data_1, data_2) {
  data_1_ex <-
    data_1[rep(1:nrow(data_1), each = nrow(data_2)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]
  data_2_ex <-
    data_2[rep(1:nrow(data_2), nrow(data_1)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]
  expanded_data <-
    data_1_ex[data_2_ex, on = 'rowid'] %>%
    .[, rowid := NULL]
  if ('tbl' %in% class(data_1)) {
    expanded_data <- as_tibble(expanded_data)
  }
  if ('rowwise_df' %in% class(data_1)) {
    expanded_data <- rowwise(expanded_data)
  }
  return(expanded_data)
}
### =======================
### Monte Carlo Simulation (for both heterogeneous field and homogeneous fields) 
### =======================
sim <- function(i_sim, field, field_raw, field_col, exp_design, coef_data, correlation_rho) {
  
  result_all <- data.table()
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
  set.seed(847259)
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
    .[, machineyerror := rnorm(nrow(data), 0, sqrt(conv_factor * var(yield_sp)))] %>%
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
    opt_N_gam <- data.table(N=seq(min(N_levels_homo),max(N_levels_homo),length=1000)) %>%
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
      .[, pi_opt := pc * yield_opt - pn * opt_N - fixed_c] %>%
      .[type=='Analysis',]
    
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
### =======================
### Generate the field ids - plot id, subplot id
### =======================
gen_field_ids <- function(field_config, field_sf) {
  
  field.row <- field_config$field.row
  field.col <- field_config$field.col
  
  # ===pre-set parameters in field_config_table
  cell_columns_in_plot <- field_config$cell_columns_in_plot
  # === number of analysis points within a plot ===#
  num_of_analysis_points_plot <- cell_columns_in_plot - buffer_in_plot
  
  temp_data <-
    CJ(row = 1:field.row, col = 1:field.col) %>%
    setorder(row, col) %>%
    # === cell_id ===#
    .[, cell_id := (row - 1) * field.col + col] %>%
    # === plot_id ===#
    .[, plot_col_id := (col - 1) %/% cell_columns_in_plot + 1] %>%
    .[, plot_row_id := (row - 1) %/% cell_rows_in_plot + 1] %>%
    .[, cell_in_plot := paste(plot_row_id, "_", plot_col_id, sep = "")] %>%
    .[, plot_id := (plot_row_id - 1) * (field.col %/% cell_columns_in_plot) + plot_col_id] %>%
    # === Analysis_Buffer_in_plot ===#
    .[, dummy := 1] %>%
    .[, cell_in_plot_in_row := cumsum(dummy), by = .(row, plot_col_id)] %>%
    .[, type := ifelse(cell_in_plot_in_row <= num_of_analysis_points_plot, "Analysis", "Buffer")] %>%
    # === subplot_id (exclude buffer zone) ===#
    .[, subplot_col_id := (col - 1) %/% cell_columns_in_subplot + 1] %>%
    .[, subplot_row_id := (row - 1) %/% cell_rows_in_subplot + 1] %>%
    .[, subplot_id := (subplot_row_id - 1) * (field.col %/% cell_columns_in_subplot) + subplot_col_id] %>%
    # === save_in_table ===#
    .[, .(
      cell_id, plot_id, subplot_id, type
    )]
  
  field_sf_with_id <- left_join(field_sf, temp_data, by = "cell_id")
  
  return(field_sf_with_id)
  
}
#' Check if cell ids and plot ids are correctly assigned
# filter(field_sf_with_id, plotid %in% 1:10) %>%
#   ggplot() +
#     geom_sf(aes(fill = cell_id)) +
#     scale_fill_viridis_c()
### =======================
### N design not using this anymore
### =======================
# fn_N_design <- function(N_levels, block_num, design){
#   #need to consider when plot_columns_in_block <>plot_rows_in_block
#   if((design=="Latin Square") & (plot_columns_in_block == plot_rows_in_block)){
#     # (1) Latin Square
#     M.latin <- rlatin(n = 1, size = num_treatment)
#     N_design <- data.table(
#       blockid=rep(1:block_num,each=num_treatment^2),
#       plot_in_block_id=rep(1:num_treatment^2,block_num),
#       N=N_levels[rep(c(t(M.latin)), times=block_num)])
#   }else{
#     M <- matrix(c(1:num_treatment),nrow=plot_rows_in_block,ncol=plot_columns_in_block,byrow=TRUE)
#     N_design <- data.table(
#       blockid=rep(1:block_num,each=plot_columns_in_block*plot_rows_in_block),
#       plot_in_block_id=rep(1:(plot_columns_in_block*plot_rows_in_block),times = block_num),
#       N=N_levels[rep(c(t(M)), times=block_num)])
#   } 	
  
#   return(N_design)
# }
#/*=================================================*/
#' # Get trial design - plot id -  N id correspondence
#/*=================================================*/
get_design_data <- function(block_col_num, block_row_num, treatment_num) {
  
  num_block_col <- floor(block_col_num)
  num_block_col_remainder <- block_col_num - floor(block_col_num)
  
  if (num_block_col >= 1) { # at least one block fits
  
    Mat_latin <- rlatin(n = 1, size = treatment_num)
    dim_mat <- matrix(1, block_row_num, num_block_col)
    base_mat <- kronecker(dim_mat, Mat_latin)
    
    #=== if there is a remainder ===#
    if (num_block_col_remainder != 0) {
      
      reduced_Mat_latin <- Mat_latin[, 1:(num_block_col_remainder * treatment_num)]
      reduced_dim_mat <- matrix(1, block_row_num, 1)
      reduced_base_mat <- kronecker(reduced_dim_mat, reduced_Mat_latin)
      whole_design <- cbind(base_mat, reduced_base_mat)
      design_plot_data <-
        data.table(n_id = c(t(whole_design))) %>%
        .[, plot_id := 1:nrow(.)]
      
    } else{
      
      whole_design <- base_mat
      design_plot_data <-
        data.table(n_id = c(t(whole_design))) %>%
        .[, plot_id := 1:nrow(.)]
      
    }
    
  } else { # only a fraction of a block fits
    
    Mat_latin <- rlatin(n = 1, size = treatment_num)
    reduced_Mat_latin <- Mat_latin[, 1:(num_block_col_remainder * treatment_num)]
    reduced_dim_mat <- matrix(1, block_row_num, 1)
    whole_design <- kronecker(reduced_dim_mat, reduced_Mat_latin)
    design_plot_data <-
      data.table(n_id = c(t(whole_design))) %>%
      .[, plot_id := 1:nrow(.)]
    
  }
  
  return(list(design_plot_data = design_plot_data, design_mat = whole_design))
}

get_strip_design_data <- function(plot_width, field_width, treatment_num){
  
  plot_num <- field_width/plot_width
  whole_design <- matrix(sample(1:treatment_num), plot_num, 1)
  design_plot_data <-   
    data.table(n_id = c(t(whole_design))) %>%
    .[, plot_id := 1:nrow(.)]
  
  return(list(design_plot_data = design_plot_data, design_mat = whole_design))
  
}
### =======================
### function to clean yield data
### =======================
sd_factor <- 4
clean_yield <- function(data, var_name) {
  temp_data <- data.table(data) %>%
    setnames(var_name, "var") %>%
    subset(var>=0)
  
  var_sd <- temp_data[
    var >= quantile(var, prob = 0.05, na.rm = TRUE) &
      var <= quantile(var, prob = 0.95, na.rm = TRUE),
    .(median = median(var, na.rm = TRUE), sd = sd(var, na.rm = TRUE))
  ]
  
  temp_data <- temp_data %>%
    .[, flag_bad := 0] %>%
    .[
      var < var_sd$median - sd_factor * var_sd$sd,
      flag_bad := 1
    ] %>%
    .[
      var > var_sd$median + sd_factor * var_sd$sd,
      flag_bad := 1
    ] %>%
    setnames("var", var_name) %>% 
    filter(flag_bad == 0)

  return(temp_data)
}
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### GWR
#/*~~~~~~~~~~~~~~~~~~~~~~*/
GWR <- function(reg_data){

  reg_data_sp <- 
    reg_data %>% 
    st_as_sf(., coords = c("X", "Y"))%>%
    as('Spatial')
  
  reg_formula <- formula(yield~N+N2)
  
  #=== gwr estimation with bw ===#
  gwr_est <- 
    gwr.basic(
      reg_formula,
      data=reg_data_sp,
      bw=100,
      kernel="gaussian",
      adaptive=T
    )
  
  # estimated beta
  gwr_beta <- 
    data.table(
      subplot_id=reg_data$subplot_id,
      b0_hat = gwr_est$SDF$Intercept,
      b1_hat = gwr_est$SDF$N,
      b2_hat = gwr_est$SDF$N2
    ) 
  # we want to calculate the opt_n_gwr for different prices with beta, outside of the GWR function
  # .[,opt_N_gwr := (b1_hat-pN/pCorn)/(-2*b2_hat)] %>%
  # #=== limit the range of opt_N_gwr ===#
  # .[,out_of_bounds:=ifelse(opt_N_gwr > max(N_levels) | opt_N_gwr < min(N_levels),1,0)] %>%
  # .[opt_N_gwr < min(N_levels), opt_N_gwr:=min(N_levels)] %>%
  # .[opt_N_gwr > max(N_levels), opt_N_gwr:=max(N_levels)]
  # #=== spatial interpolate the buffer zone cells ===#
  # {	
  # 	# cell level data
  # 	sdata <- data[, .(cell_id, subplotid)] %>% 
  # 		gwr_beta[.,on='subplotid'] %>% 
  # 		left_join(field[, c('cell_id','type')], ., by='cell_id') %>%
  # 		cbind(., st_coordinates(st_centroid(.))) %>%
  # 		data.table()
  
  # 	# effective sample points
  # 	PT <- sdata[type=='Analysis', .(type, X, Y, cell_id, subplotid, b0_hat, b1_hat, b2_hat)] # took out opt_N_gwr
  # 	coordinates(PT) = ~X+Y      # spatial points
  # 	# whole field grids
  # 	PX <- sdata[, .(type, X, Y, cell_id, subplotid)]
  # 	coordinates(PX) = ~X+Y  	# spatial points
  # 	gridded(PX) = TRUE      	# spatial pixel
  
  # 	# IDW interpolated opt N rates
  # 	PX$b0_hat = gstat::idw(b0_hat~1, PT, PX)$var1.pred
  # 	PX$b1_hat = gstat::idw(b1_hat~1, PT, PX)$var1.pred
  # 	PX$b2_hat = gstat::idw(b2_hat~1, PT, PX)$var1.pred
  # 	# -> the in-sample points are not interpolated, which
  # 	#	is great, but not sure why
  
  # 	# aggregate N rates back to subplot level
  # 	gwr_beta_full <- PX@data %>% data.table() %>%
  # 		.[, .(b0_hat=mean(b0_hat), b1_hat=mean(b1_hat), b2_hat=mean(b2_hat)), by=subplotid]
  
  # 	# update gwr_beta
  
  # 	gwr_beta <- if(all(PT$b0_hat==PX$b0_hat[PX$type=='Analysis'])){gwr_beta_full} 
  # 	else {gwr_beta <- gwr_beta}
  # }
  return(gwr_beta)
}
# left_join(field, data[,.(cid,opt_N)],by='cid') %>%
#   select(opt_N) %>%
#   plot()
