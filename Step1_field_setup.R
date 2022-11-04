## Field Configuration
## Load packages
library("raster")
library("dplyr")
library("sf")
library("data.table")
library("magrittr")
library("magic")
library("here")
source(here("./Code/Codes/Functions/functions_for_simulation.R"))
# /*=================================================*/
##                      Parameters
# /*=================================================*/

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~
# === plot length ===#
# I am setting the cell numbers so that we would have integer number of plots in one row of the field
plot_len_ls <- c(6 + 6, 6 + 18, 6 + 42, 999) # 1.79m per cell # currently 5 plot length
# === field size ===#
cell_columns_in_field_ls <- c(144, 288, 576) # three different field sizes
cell_rows_in_field_ls <- c(288)
cell_field <-
  CJ(
    cell_columns_in_field = cell_columns_in_field_ls,
    cell_rows_in_field = cell_rows_in_field_ls
  )
saveRDS(cell_field, here("Data", "cell_field.rds"))
# === number of treatment levels ===#
num_treatment_ls <- c(4, 6, 8) # these number of treatment rates alligns well with the field dimension
# === yield monitor accuracy ===#
correlation_rho_ls <- c(0.7,0.9)
# /*=================================================*/
##                      Field layout
# /*=================================================*/
#' cell: basic unit
#' plot: N application unit

# ~~~~~~~~~~~~~~~
# sizes of field
# ~~~~~~~~~~~~~~~
# === in cells ===#
cell_rows_in_plot <- 12 # depends on yield monitor
buffer_in_plot <- 6
# === unit of analysis ===# ("subplot")(same across any design)
cell_rows_in_subplot <- 12
cell_columns_in_subplot <- 6
# === in meters ===#
cell <- 1.79

#--------------------------
# Experimental settings
#--------------------------
#' field size
#' plot length
#' number of treatments
# === field set up summary ===#
exp_designs <- -data.table(
  "field.size" = numeric(),
  "field.cell" = numeric(),
  "field.col" = numeric(),
  "field.lgth" = numeric(),
  "field.row" = numeric(),
  "field.wdth" = numeric(),
  "correlation_rho" = numeric(),
  "cell_columns_in_plot" = numeric(),
  "num_treatment" = numeric(),
  "plot_columns_in_field" = numeric(),
  "plot_rows_in_field" = numeric()
)

for (i in 1:nrow(cell_field)) {
  cell_columns_in_field <- cell_field[i, cell_columns_in_field] # field size
  cell_rows_in_field <- cell_field[i, cell_rows_in_field]

  for (j in 1:length(plot_len_ls)) {
    cell_columns_in_plot <- ifelse(plot_len_ls[j] == 999, cell_columns_in_field, plot_len_ls[j]) # plot length
    plot_columns_in_field <- cell_field[i, cell_columns_in_field] %/% cell_columns_in_plot # field size
    plot_rows_in_field <- cell_field[i, cell_rows_in_field] %/% cell_rows_in_plot
    field.col <- cell_columns_in_field # field columns in cells
    field.row <- cell_rows_in_field # field rows in cells
    field.cell <- field.col * field.row
    field.size <- round(1.79^2 * field.cell * 10^(-4),2)
    # basic cell (meters)
    field.lgth <- field.col * cell
    field.wdth <- field.row * cell # field

    for (k in 1:length(num_treatment_ls)) {
      num_treatment <- num_treatment_ls[k]
        for (c in 1:length(correlation_rho_ls)){
          correlation_rho <- correlation_rho_ls[c]
          field_table_new <- data.table(
            field.size, field.cell, field.col, field.lgth,
            field.row, field.wdth, correlation_rho, cell_columns_in_plot, num_treatment,
            plot_columns_in_field, plot_rows_in_field
          )
          exp_designs <- rbind(exp_designs, field_table_new)
        }
    }
  }
}

saveRDS(exp_designs, here("Data", "exp_designs.rds"))
# exp_designs <- readRDS(here("Data","exp_designs.rds"))
#/*=================================================*/
#' # Create fields of different size
#/*=================================================*/
for (f in 1:nrow(cell_field)) {

  field.row <- cell_field[f][[2]]
  field.col <- cell_field[f][[1]]
  field.cell <- field.row * field.col
  field.lgth <- field.col * cell
  field.wdth <- field.row * cell

  Mf <- 
    matrix(
      1:field.cell, 
      nrow = field.row, 
      ncol = field.col, 
      byrow = TRUE
    )

  # === matrix -> raster
  field_raster <- raster(Mf)
  extent(field_raster) <- extent(0, field.lgth, 0, field.wdth)
  names(field_raster) <- "cell_id"

  # === rater -> polygon sf
  field_sf <- 
    as(field_raster, "SpatialPolygonsDataFrame") %>%
    st_as_sf() %>% 
    cbind(., st_coordinates(st_centroid(.)))   

  saveRDS(field_sf, paste0("Data/field_col_", field.col, ".rds"))

}
#/*=================================================*/
#' # Define the treatment matrix for the different settings
#/*=================================================*/
# exp_designs <- readRDS(here("Data", "exp_designs.rds"))
# library(gridExtra)
# pdf("data_output.pdf", height=11, width=15)
# grid.table(exp_designs)
# dev.off()

set.seed(893253)

# head(temp$design_plot_data, 50)
# temp <-   
#   get_design_data(
#     block_col_num = trial_designs[4, ]$block_col_num, 
#     block_row_num = trial_designs[4, ]$block_row_num, 
#     treatment_num = trial_designs[4, ]$num_treatment
#   )

trial_designs <- 
  exp_designs %>%
  .[, block_col_num := field.col / (cell_columns_in_plot * num_treatment)] %>%
  .[, block_row_num := field.row / (cell_rows_in_plot * num_treatment)] %>% 
  rowwise() %>% 
  mutate(trial_design = list(
    if (field.col != cell_columns_in_plot) {
      #=== non-strip trial ===#
      get_design_data(block_col_num, block_row_num, num_treatment)$design_plot_data
      # get_design_data(block_col_num, block_row_num, num_treatment)$design_mat
    } else {
      #=== strip trial ===#
      get_strip_design_data(cell_rows_in_plot, field.row, num_treatment)$design_plot_data
      # get_strip_design_data(cell_rows_in_plot, field.row, num_treatment)$design_mat
    }  
  ))

saveRDS(trial_designs, "Data/trial_design.rds")

# head(trial_designs$trial_design[[4]], 50)