#===================================
# Preparation
#===================================
#--- loading packages ---#
library("data.table")
library("dplyr")
library("sf")
library("gstat")
library("magrittr")
library("here")
library("ggplot2")
source(here("./Functions/functions_for_simulation.R"))

#--------------------------
# set number of simulations
#--------------------------
num_sim_p <- 500
gstat_model <- "Sph"
sp_range <- 600
# pCorn <- 4.1925 / 25.40 # $/kg #($4.1925 per bushel, 1 bushel = 25.40kg)
# pN <- 0.47 / 0.453592 # per kg (#0.47 per lb n /0.453592)

cell_field <- readRDS(here("Data", "cell_field.rds"))

for (i in 1: nrow(cell_field)){

  temp_field <- cell_field[i, ]
  
  field <- readRDS(paste0("Data/field_col_", temp_field$cell_columns_in_field, ".rds"))
      
  xy <- 
    st_centroid(field) %>%
    st_coordinates() %>%
    data.table() %>%
    .[, cell_id := field$cell_id]

  source(here("./Functions/parameter_generation.R"))
  
  saveRDS(
    par_data,
    here(
      "Data", 
      paste0(
        'field_col_', 
        temp_field$cell_columns_in_field,
        '_coef.rds'
      )
    )
  )

}

