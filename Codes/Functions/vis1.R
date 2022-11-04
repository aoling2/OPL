
all_hetero_result <- data.table()
for (i in 1:length(sim_hetero_results_ls)) {
  all_hetero_result <- rbind(all_hetero_result, sim_hetero_results_ls[i][[1]])
}

all_sim_result <- readRDS(here("./Results/result.rds"))

pi_data <- all_sim_result %>%
  .[order(field_size, plot_length, num_treatment, yield_accuracy),] %>%
  print(digits=2)

mean_data <- pi_data %>%
  #===take average across simulations
  .[, .(pi_opt=mean(pi_opt),
        pi_gwr=mean(pi_gwr),
        pi_brf=mean(pi_brf),
        opt_N=mean(opt_N),
        opt_N_gwr=mean(opt_N_gwr),
        opt_N_rf_perfect=mean(opt_N_rf_perfect),
        pi_opt_homo=mean(pi_opt_homo),
        pi_opt_gam=mean(pi_opt_gam),
        opt_N_homo=mean(opt_N_homo),
        opt_N_gam=mean(opt_N_gam)), by=c("field_size", "plot_length", "num_treatment", "yield_accuracy")] %>%
  print()
