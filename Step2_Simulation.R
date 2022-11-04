# ===================================
# Preparation
# ===================================
#--- loading packages ---#
library("data.table")
library("dplyr")
library("parallel")
library("ggplot2")
library("sf")
library("ggplot2")
library("gstat")
library("magrittr")
library("tictoc")
library("GWmodel")
library("here")
library("grf")
library("Rcpp")
library("mgcv")
source(here("./Code/Codes/Functions/functions_for_simulation.R"))
# Start the clock!

tic()

# /*=================================================*/
#' # Set up
# /*=================================================*/
# === number of simulations ===#
num_sim <- 500

# === load in all experimental design settings ===#
trial_designs <- readRDS(here("Code/Data/trial_design.rds"))
exp_designs <- readRDS(here("Code/Data/exp_designs.rds"))
# === unit size ===# 
buffer_in_plot <- 6
cell_rows_in_plot <- 12 
cell_rows_in_subplot <- 12
cell_columns_in_subplot <- 6
# === spatial feature ===# 
gstat_model <- "Sph"
sp_range <- 600
#--------------------------
# parameters unchanged for all simulations
#--------------------------
pCorn <- 4.1925 / 25.40 # $/kg #($4.1925 per bushel, 1 bushel = 25.40kg)
pN <- 0.47 / 0.453592 # per kg (#0.47 per lb n /0.453592)

# === fixed cost ===#
fixed_c <- c(1112) # $/ha $450/acre Taro uses $550/acre

#--- define the range of Pc and Pn ---#
# $ per bushel corn ranges from 1 to 5.5 (/25.40)=0.04 to 0.22
# $ per lb nitrogen ranges from 0.05 to 0.48 (/0.453592)=0.10 to 1
# 6.277592 is the current Pn/Pc

Pc_ls <- c(pCorn, 1, 1)
Pn_ls <- c(pN,3, 10)

price_eval <- cbind(Pc_ls, Pn_ls) %>%
  data.table() %>%
  rename(Pc = 'Pc_ls', Pn = 'Pn_ls') %>%
  mutate(price_ratio := Pn / Pc)%>%
  round(.,2)

price_eval <- price_eval[1,]

#-------------------------#
#  Nitrogen range (for homo fields)
#-------------------------#
N_homo_4 <- c(80, 160, 200, 240)*1.12085
N_homo_6 <- c(90, 120, 150, 180, 210, 240)*1.12085
N_homo_8 <- c(80, 100, 120, 160, 180, 200, 220, 240)*1.12085
#/*=================================================*/
#' # Simulations
#/*=================================================*/
sim_results_ls <- list()
for (i in 1:24) {
{
  temp_design <- trial_designs[i, ]
  field_col <- temp_design$field.col

  field_raw <- readRDS(paste0("Code/Data/field_col_", field_col, ".rds"))

  coef_data <- 
    readRDS(
      here(
        "Coode/Data", 
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

  # /*----------------------------------*/
  #' ## Simulations
  # /*----------------------------------*/
  # === simulations ===#
}
print(i)
{
    sim_results <-
      mclapply(1:num_sim,
               function(x) {
                 sim(
                   i_sim = x,
                   field = field,
                   field_raw = field_raw,
                   field_col = field_col,
                   exp_design = temp_design,
                   coef_data = coef_data,
                   correlation_rho = correlation_rho
                 )
               },
               mc.cores = 8
      ) %>%
      rbindlist() %>%
      .[, field_size := temp_design$field.size] %>%
      .[, yield_accuracy := yield_accuracy] %>%
      .[, plot_length := temp_design$cell_columns_in_plot] %>%
      .[, num_treatment := temp_design$num_treatment] %>%
      .[, i_sim := i_sim] %>%
      .[, pC_pN := paste(pCorn, "_", pN, sep = "")] %>%
      .[, pCorn := pCorn] %>%
      .[, pN := pN] %>%
      .[, opt_N_gwr := opt_N_gwr] %>%
      .[, pi_gwr := pi_gwr] %>%
      .[, opt_N_rf_perfect := opt_N_rf_perfect] %>%
      .[, pi_brf := pi_brf] %>%
      .[, opt_N := opt_N] %>%
      .[, pi_opt := pi_opt] %>%
      .[, opt_N_gam := opt_N_gam] %>%
      .[, pi_opt_gam := pi_opt_gam] %>%
      .[, opt_N_homo := opt_N_homo] %>%
      .[, pi_opt_homo := pi_opt_homo]

    sim_results_ls[[i]] <- sim_results
  }
}

# Save the results
all_sim_result <- data.table()
for (i in 1:length(sim_results_ls)) {
  all_sim_result <- rbind(all_sim_result, sim_results_ls[i][[1]])
}

saveRDS(all_sim_result, here("./Code/Results/result_500_1_24.rds"))
# Stop the clock
toc()
system("say Aolin, just finished!")

