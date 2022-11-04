i <- 24
for (i in 1:nrow(exp_designs)){
  field.col <- exp_designs[i, ]$field.col
  field.row <- exp_designs[i, ]$field.row
  cell_columns_in_plot <- exp_designs[i, ]$cell_columns_in_plot
  cell_rows_in_plot <- 12
  num_treatment <- exp_designs[i, ]$num_treatment
  block_col_num <- field.col / (cell_columns_in_plot * num_treatment)
  block_row_num <- field.row / (cell_rows_in_plot * num_treatment)


  get_design_data(block_col_num, block_row_num, num_treatment)$design_plot_data
  get_design_data(block_col_num, block_row_num, num_treatment)$design_mat
  #=== strip trial ===#
  get_strip_design_data(cell_rows_in_plot, field.row, num_treatment)$design_plot_data
  get_strip_design_data(cell_rows_in_plot, field.row, num_treatment)$design_mat
}
get_design_data <- function(block_col_num, block_row_num, treatment_num) {
  treatment_num <- num_treatment
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
}
  
trial_designs <- readRDS(here("Data/trial_design.rds"))
test <- trial_designs[24,]$trial_design[[1]]
