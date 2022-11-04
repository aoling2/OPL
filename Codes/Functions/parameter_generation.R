{
# y = -b2 * x2 + b1 * x - b0
# === N_star ===#
N_star <- gen_coefs(xy, mean = 200, psill = 1000, range = sp_range, coef_name = "N_star", nsim = num_sim_p) %>%
  # >>> normalize <<<
  .[,sd_b:=sd(N_star),by=sim] %>%
  .[,mean_b:=mean(N_star),by=sim] %>%
  .[,p:=pnorm(N_star,mean=mean_b, sd=sd_b)] %>%
  .[,N_star:=100+p*150] %>%
  .[, c("cell_id", "sim", "N_star")]
# N_star%>%
#   .[,N_star]%>%
#   hist()
# === ymax ===#
ymax <- gen_coefs(xy, mean = 15000, psill = 3000000, range = sp_range, coef_name = "ymax", nsim = num_sim_p) %>%
  # >>> normalize <<<
  .[,sd_b:=sd(ymax),by=sim] %>%
  .[,mean_b:=mean(ymax),by=sim] %>%
  .[,p:=pnorm(ymax,mean=mean_b, sd=sd_b)] %>%
  .[,ymax:=8000+p*8000] %>%
  .[, c("cell_id", "sim", "ymax")]
# ymax%>%
#   .[,ymax]%>%
#   hist()
# === b0 ===#
b0 <- gen_coefs(xy, mean = 6000, psill = 200000, range = sp_range, coef_name = "b0", nsim = num_sim_p) %>%
  .[ymax, on = c("sim", "cell_id")] %>%
  # >>> normalize <<< (b0 needs to be < ymax)
  .[,sd_b:=sd(b0),by=sim] %>%
  .[,mean_b:=mean(b0),by=sim] %>%
  .[,p:=pnorm(b0,mean=mean_b, sd=sd_b)] %>%
  .[,b0:=3000+p*4000] %>%
  .[, c("cell_id", "sim", "b0")]
# b0%>%
#   .[,b0]%>%
#   hist()
# === b1, b2 ===#
coef_data <- N_star %>%
  .[ymax, on = c("sim", "cell_id")] %>%
  .[b0, on = c("sim", "cell_id")] %>%
  # === derive b1, b2 from b0, ymax, and N_star
  .[, b2 := (ymax - b0) / N_star^2] %>%
  .[, b1 := 2 * N_star * b2] 
  # .[, opt_N := (- pN / pCorn + b1) / (2 * b2)] 
# 
# coef_data%>%
#   .[,b1]%>%
#   hist()
# coef_data%>%
#   .[,b0]%>%
#   hist()
# coef_data%>%
#   .[,b2]%>%
#   hist()
# 
# coef_data%>%
#   .[,opt_N]%>%
#   hist()
#b1 15 to 150
# yield_function <- function(x, b0, b1, b2){
#   y <- -b2 * x^2 + b1 * x + b0
# }
# 
# data_test <- 
#   coef_data %>%
#   .[sim == 1, ] %>%
#   .[, c("cell_id", 'b0', 'b2', 'b1')] %>%
#   expand_grid_df(., data.table(x = seq(0, 300, by = 10))) %>%
#   .[, yield := yield_function(x, b0, b1, b2)]
# 
# ggplot(data_test[b0 < 4000, ][1:10000, ]) +
#   geom_line(aes(y = yield, x = x, group = cell_id)) +
#   xlim(0, 200) +
#   ylim(0, NA)
# 
# ymax = b0 + b1N_star + b2N_star2
# -b0-b2N_star^2 = ymax - b1N_star
# b1 = (ymax-b0-b2N_star^2)/N_star
# b2 = (ymax-b0-b1N_star)/N_star2

  # y = (b0 - alphaN_star^2) - alpha x^2 + 2alphaN_star x
  # b1 = 2alpha N_star
  # b2 = -alpha
  # N_star = - b1/(2b2)

  # ymax = (b0 - alphaN_star^2) - alpha N_star^2 + 2alphaN_star^2
# === plateau/Nk ===# (simply use N_star as Nk)
coef_data <- coef_data %>%
  .[, Nk := N_star] %>%
  .[, plateau := ymax]%>%
  .[, c("cell_id", "sim", "b0", "b1", "b2", "Nk", "plateau")]

m_error <- gen_errors(mean = 0, psill = 0.015, range = sp_range, coef_name = "m_error", nsim = num_sim_p) %>%
  # >>> normalize <<<
  .[, min_b := min(m_error), by = sim] %>%
  .[, max_b := max(m_error), by = sim] %>%
  .[, p := punif(m_error, min_b, max_b)] %>%
  .[, m_error := -0.3 + p * 0.6] %>%
  .[, c("cell_id", "sim", "m_error")]

# ====cell-level N errors====#
N_error <- gen_errors(mean = 0, psill = 0.2, range = 50, coef_name = "N_error", nsim = num_sim_p) %>%
  # >>> normalize <<<
  .[, min_b := min(N_error), by = sim] %>%
  .[, max_b := max(N_error), by = sim] %>%
  .[, p := punif(N_error, min_b, max_b)] %>%
  .[, N_error := -1 + p * 2] %>%
  .[, c("cell_id", "sim", "N_error")]

# === combine all ===#
par_data <- m_error[coef_data, on = c("cell_id", "sim")] %>% N_error[., on = c("cell_id", "sim")]
}