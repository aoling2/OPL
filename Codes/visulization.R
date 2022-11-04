## ----echo = F-----------------------------------------------------------------
library(knitr)
library(here)
library(officedown)
library(officer)

opts_chunk$set(
  fig.align = "center",
  fig.retina = 5,
  warning = F,
  message = F,
  cache = T,
  echo = F,
  fig.cap = TRUE
)


## ----echo = F-----------------------------------------------------------------
#=== Packages ===#
library(sf)
library(raster)
library(ggplot2)
library(data.table)
library(magrittr)
library(viridis)
library(dplyr)
library(data.table)
library(modelsummary)
library(patchwork)
library(flextable)
library(here)
library(gstat)
library(latex2exp)
source(here("./Code/Codes/Functions/functions_for_simulation.R"))


## ----fig.id = "Nstar-map", cache = TRUE, fig.cap = "Simulated Spatial Distribution of true EONR (Variogram Range=400 meters)"----
# raster map: plot the map when you first generate the data using gstat() function
num_sim_p <- 1
gstat_model <- "Sph"
sp_range <- 400
field <- readRDS(here(paste0("Code/Data/field_col_", 144, ".rds")))
xy <- 
  st_centroid(field) %>%
  st_coordinates() %>%
  data.table() %>%
  .[, cell_id := field$cell_id]
cell_data <- gen_coefs(xy, mean=200,psill=1000,range=sp_range,coef_name='N_star',nsim=num_sim_p) %>%
    #>>> normalize <<<
    .[,sd_b:=sd(N_star),by=sim] %>%
    .[,mean_b:=mean(N_star),by=sim] %>%
    .[,p:=pnorm(N_star,mean=mean_b, sd=sd_b)] %>%
    .[,N_star:=100+p*150] %>%
    .[,c('cell_id','sim','X','Y','N_star')]
# transfer X Y
cell_data <- cell_data %>%
  .[, X := 0.895 + as.numeric(as.factor(X)) * 1.79 - 1.79] %>%
  .[, Y := 0.895 + as.numeric(as.factor(Y)) * 1.79 - 1.79]

N_star_raster <- ggplot() +
                geom_raster(data = cell_data[sim==1, c("X","Y","N_star")] , 
                            aes(x = X, y = Y, fill = N_star)) +
                scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"))
                # scale_fill_gradientn(colours=c("chocolate1","chocolate2","sienna2","sienna1","palegreen2","palegreen3"))

N_star_raster


## ----fig.id = "yield-response", cache = TRUE, fig.cap = "The quadratic plateau of the estimated yield response to nitrogen with the estimated parameters"----
a = 216.25625858
b = 10.10459998
c = -0.06728008

low <- function(x){
  (a + b * x + c * I(x^2)) * (x <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c)
}

plot_data <- data.table(x = seq(8, 82, length = 1000)) %>% 
 .[,y_VIlow := low(x)] %>% 
  melt(id.var = "x") %>% 
  .[, type := case_when(
    variable == "y_VIlow" ~ "VI_low"
  )]

find_point_x <- plot_data%>%filter(variable == 'y_VIlow', x >= 62, x<=62.5)

df_conceptual <- data.frame(
  xl = find_point_x[6, x],
  yl = find_point_x[6, value]
)

slope_conceptual = 0.15

df_conceptual <- df_conceptual%>%
  data.table()%>%
  .[, `:=` (
    xl1 = xl + 11,
    xl2 = xl - 14,
    yl1 = yl + slope_conceptual*11*12.5+1,
    yl2 = yl - slope_conceptual*14*12.5+1
  )]

df_conceptual_1 <- data.frame(
  x = df_conceptual[, c(xl1, xl2)],
  y = df_conceptual[, c(yl1, yl2)],
  group = c('VI_low', 'VI_low')
)

#--- conceptual steps figure ---#
conceptual_yield <- ggplot()+
  geom_line(aes(x = x, y = value, color = type), size = 0.9, data = plot_data, se = FALSE)+
  scale_color_manual(
    values = c(
      "VI_low" = "#2E86C1"
    )
  )+
  geom_linerange(aes(x=8, y=NULL, ymin=240, ymax=740))+
  geom_linerange(aes(x=NULL, y=240, xmin=8, xmax=85))+
  labs(x = 'Side-dressing N Rate', y ='Yield') +
  geom_point(aes(x = xl, y = yl), data = df_conceptual)+
  #--- tangent line ---#
  geom_line(aes(x = x, y = y, group = group), color = '#000000', size = 0.4, data = df_conceptual_1)+
  #--- angle h-line ---#
  geom_segment(aes(x = df_conceptual$xl2, y = df_conceptual$yl2, xend = df_conceptual$xl, yend = df_conceptual$yl2), size = 0.4, linetype = 'solid')+
  #--- EONR vline ---#
  geom_segment(aes(x = df_conceptual$xl, xend = df_conceptual$xl, y = df_conceptual$yl, yend = 240), linetype = 'dashed', size = 0.3)+
  #--- curve ---#
  geom_curve(aes(x = df_conceptual$xl - 5.8, xend = df_conceptual$xl - 5.5, y = df_conceptual$yl - 10, yend = df_conceptual$yl2), curvature = -0.3)+
  #--- text for the ratio ---#
  geom_text(aes(x = df_conceptual$xl2, y =510, label = TeX("$\\frac{P_N}{P_C}$", output = "character")), parse = TRUE, stat = 'unique', size = 5)+
  geom_segment(aes(x = df_conceptual$xl2+2, xend = df_conceptual$xl - 6, y = df_conceptual$yl2- 35, yend = df_conceptual$yl - 18), 
               arrow = arrow(length = unit(0.18, 'cm')), size = 0.3)+
  #--- text for estimated yield response curve at site1 ---#
  annotate('text', x = 78, y =500, label = 'estimated yield response curve', size = 5, color = '#2E86C1', family = 'Times')+
  annotate('text', x = 85, y = 595, parse = TRUE, label = (TeX(r'($$\textit{\tilde{f}(N^{i}})$$)')), color = "#2E86C1", size = 5.4, family = "Times", fontface = 5)+
  geom_segment(aes(x = 73, xend = 78, y = 510, yend = 590), arrow = arrow(length = unit(0.15, 'cm')), size = 0.3, color = '#2E86C1')+
  annotate('text', x = df_conceptual$xl, y = 240-19, label = (TeX(r'($$\hat{N}^{i*}$$)')), size = 5, color = '#2E86C1', family = 'Times')+
  annotate('text', x = 81, y = 240-19, label = 'N application rate', size = 5.5, family = 'Times')+
  annotate('text', x = 8, y = 750, label = 'Yield', size = 5.5, family = 'Times')+
  coord_fixed(0.1)+
  xlim(8, 92)+
  theme_void()+
  theme(legend.position = "none")
conceptual_yield


## ----fig.id = "treatment-map", cache = TRUE, fig.cap = "The treatment map of the 2017 checkerboard trial"----
treatment <- read_sf(here('./Code/Data/trial-design.shp'))
treatment_map <- ggplot(treatment) +
  geom_sf(aes(fill = (NRATE_Gal3*11.06*0.32+192)*1.12085)) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme_bw() +
  theme_void() +
  labs(fill = "N rates (kg/ha)") +
  theme(legend.position = "bottom") +
  theme(
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    legend.key.width = unit(0.6, "cm")
  )

treatment_map


## ----load data, include=FALSE-------------------------------------------------
all_sim_result <- readRDS(here('./Code/Results/result_500.rds'))

pi_data <- all_sim_result %>%
  .[order(field_size, 
          plot_length, 
          num_treatment,
          yield_accuracy),] %>%
  print(digits=2)

mean_data <- pi_data %>%
  #===take average across simulations
  .[, .(pi_opt=mean(pi_opt),
        pi_gwr=mean(pi_gwr),
        pi_brf=mean(pi_brf),
        pi_gwr_dif=mean(pi_gwr-pi_opt),
        pi_brf_dif=mean(pi_brf-pi_opt),
        pi_gwr_sd=sd(pi_gwr-pi_opt),
        pi_brf_sd=sd(pi_brf-pi_opt),
        u_gwr=2.5*mean((pi_gwr)^(0.4)),
        u_brf=2.5*mean((pi_brf)^(0.4)),
        opt_N=mean(opt_N),
        opt_N_gwr=mean(opt_N_gwr),
        sd_N_gwr=sd(opt_N_gwr-mean(opt_N)),
        opt_N_rf_perfect=mean(opt_N_rf_perfect),
        sd_N_brf=sd(opt_N_rf_perfect-mean(opt_N)),
        pi_opt_homo=mean(pi_opt_homo),
        pi_opt_gam=mean(pi_opt_gam),
        pi_gam_dif=mean(pi_opt_gam-pi_opt_homo),
        pi_gam_sd=sd(pi_opt_gam-pi_opt_homo),
        u_gam=2.5*mean((pi_opt_gam)^(0.4)),
        opt_N_homo=mean(opt_N_homo),
        opt_N_gam=mean(opt_N_gam), 
        sd_N_gam=sd(opt_N_gam-mean(opt_N_homo))),
    by=c("field_size", 
         "yield_accuracy",
         "plot_length", 
         "num_treatment"
         )] %>%
  print(digits=2)


## ----var names----------------------------------------------------------------
var_names <- c(
  `0.7` = "yield accuracy level = 0.7",
  `0.9` = "yield accuracy level = 0.9",
  `13.29` = "field size = 13.29",
  `26.58` = "field size = 26.58",
  `4` = "4 treatments",
  `6` = "6 treatments",
  `8` = "8 treatments",
  `12` = "plot length = 12",
  `24` = "plot length = 24",
  `48` = "plot length = 48",
  `144` = "plot length = strip",
  `288` = "plot length = strip"
)


## ----tab.id = "data-summary-0", tab.cap = "All combinations of experimental designs"----
exp_design <- readRDS(here('./Code/Data/exp_designs.rds'))
data_summary_0 <- exp_design %>%
  filter(field.size<53.15) %>%
  dplyr::select(field.size,
         field.lgth,
         cell_columns_in_plot,
         num_treatment,
         correlation_rho,) %>%
  flextable()%>%
  theme_vanilla() %>%
  set_header_labels(
    field.size = "field size",
    field.lgth = 'field length',
    cell_columns_in_plot = "plot length",
    num_treatment = "number of treatment",
    correlation_rho = "yield accuracy"
  ) %>%
  merge_v(j = c("field.size", 'field.lgth', "cell_columns_in_plot","num_treatment")) %>%
  italic(j = 1) %>%
  align_nottext_col(., align = "center") %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.4, part = "all") %>%
  autofit()

data_summary_0


## ----tab.id = "data-summary-1", tab.cap = "GWR Simulation results of different experimental designs on a heterogeneous field of 13.29 hectares"----
data_summary_1 <- mean_data %>%
  filter(., field_size == 13.29) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_gwr_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_gwr_dif = "profit deficit estimated profit by gwr"
  ) %>%
  highlight(
    i = ~ pi_gwr_dif == min(pi_gwr_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_gwr_dif == max(pi_gwr_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_gwr_dif'), digits = 2) %>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_1


## ----tab.id = "data-summary-2", tab.cap = "GWR Simulation results of different experimental designs on a heterogeneous field of 26.58 hectares"----
data_summary_2 <- mean_data %>%
  filter(., field_size == 26.58) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_gwr_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_gwr_dif = "profit deficit estimated profit by gwr"
  ) %>%
  highlight(
    i = ~ pi_gwr_dif == min(pi_gwr_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_gwr_dif == max(pi_gwr_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_gwr_dif'), digits = 2) %>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_2


## ----tab.id = "data-summary-3", tab.cap = "BRF Simulation results of different experimental designs on a heterogeneous field of 13.29 hectares"----
data_summary_3 <- mean_data %>%
  filter(., field_size == 13.29) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_brf_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_brf_dif = "profit deficit estimated profit by brf"
  ) %>%
  highlight(
    i = ~ pi_brf_dif == min(pi_brf_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_brf_dif == max(pi_brf_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_brf_dif'), digits = 2)%>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_3


## ----tab.id = "data-summary-4", tab.cap = "BRF Simulation results of different experimental designs on a heterogeneous field of 26.58 hectares"----
data_summary_4 <- mean_data %>%
  filter(., field_size == 26.58) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_brf_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_brf_dif = "profit deficit estimated profit by brf"
  ) %>%
  highlight(
    i = ~ pi_brf_dif == min(pi_brf_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_brf_dif == max(pi_brf_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_brf_dif'), digits = 2)%>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_4


## ----tab.id = "data-summary-5", tab.cap = "GAM Simulation results of different experimental designs on a homogeneous field of 13.29 hectares"----
data_summary_5 <- mean_data %>%
  filter(., field_size == 13.29) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_gam_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_gam_dif = "profit deficit estimated profit by gam"
  ) %>%
  highlight(
    i = ~ pi_gam_dif == min(pi_gam_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_gam_dif == max(pi_gam_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_gam_dif'), digits = 2)%>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_5


## ----tab.id = "data-summary-6", tab.cap = "GAM Simulation results of different experimental designs on a homogeneous field of 26.58 hectares"----
data_summary_6 <- mean_data %>%
  filter(., field_size == 26.58) %>%
  dplyr::select(plot_length,
         num_treatment,
         yield_accuracy,
         pi_gam_dif
  ) %>%
  flextable()%>%
  merge_v(j = c("plot_length","num_treatment")) %>%
  set_header_labels(
    plot_length = "plot length",
    num_treatment = 'number of treatment',
    yield_accuracy = "yield accuracy",
    pi_gam_dif = "profit deficit estimated profit by gam"
  ) %>%
  highlight(
    i = ~ pi_gam_dif == min(pi_gam_dif),
    j =c(1,2,4),
    color = "red"
  ) %>%
  highlight(
    i = ~ pi_gam_dif == max(pi_gam_dif),
    j =c(1,2,4),
    color = "green"
  ) %>%
  autofit() %>%
  border_inner_h(border = fp_border(color = "black")) %>%
  border_inner_v(border = fp_border(color = "black")) %>%
  line_spacing(space = 0.5, part = "all") %>%
  align_text_col(., align = "justify") %>%
  align_nottext_col(., align = "center") %>%
  colformat_double(j = c('pi_gam_dif'), digits = 2)%>%
  add_footer(plot_length = "Note: The design that achieved the best outcome for each variable is highlighted in green, and the worst outcome in red.") %>%
  merge_at(j = 1:4, part = "footer")

data_summary_6


## ----fig.id = "num-treatment-gwr", cache = TRUE, fig.cap = "Impact of number of treatment on profit estimated by GWR"----
g_num_treats_gwr <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_gwr_dif:=pi_gwr-pi_opt) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_gwr_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=.) +
  geom_line(aes(x=num_treatment, y=mean)) +
  geom_point(aes(x=num_treatment, y=mean)) +
  scale_x_continuous(breaks = c(4,6,8)) +
  facet_grid(plot_length~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Number of treatment') + 
  ylab('Estimated Deficit Profit \n by Geographically Weighted Regression (Kg/Ha)')

g_num_treats_gwr


## ----fig.id = "num-treatment-brf", cache = TRUE, fig.cap = "Impact of number of treatment on profit estimated by BRF"----
g_num_treats_brf <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_brf_dif:=pi_brf-pi_opt) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_brf_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=.) +
  geom_line(aes(x=num_treatment, y=mean)) +
  geom_point(aes(x=num_treatment, y=mean)) +
  scale_x_continuous(breaks = c(4,6,8)) +
  facet_grid(plot_length~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Number of treatment') + 
  ylab('Estimated Deficit Profit \n by Balanced Random Forest (Kg/Ha)')

g_num_treats_brf


## ----fig.id = "num-treatment-gam", cache = TRUE, fig.cap = "Impact of number of treatment on profit estimated by GAM"----
g_num_treats_gam <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_gam_dif:=pi_opt_gam-pi_opt_homo) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_gam_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=.) +
  geom_line(aes(x=num_treatment, y=mean)) +
  geom_point(aes(x=num_treatment, y=mean)) +
  scale_x_continuous(breaks = c(4,6,8)) +
  facet_grid(plot_length~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Number of treatment') + 
  ylab('Estimated Deficit Profit \n by Generalized Additive Model (Kg/Ha)')

g_num_treats_gam


## ----fig.id = "plot-length-gwr", cache = TRUE, fig.cap = "Impact of plot length on profit estimated by GWR"----
g_plot_length_gwr <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_gwr_dif:=pi_gwr-pi_opt,
         plot_length:=factor(plot_length,levels=c("12","24","48","Strip"))) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_gwr_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=., aes(x=plot_length, y=mean, group=num_treatment)) +
  geom_line() +
  geom_point() +
  scale_x_discrete("plot_length", labels=c("12","24","48","Strip")) +
  facet_grid(num_treatment~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Plot Length') + 
  ylab('Estimated Deficit Profit \n by Geographically Weighted Regression (Kg/Ha)')

g_plot_length_gwr


## ----fig.id = "plot-length-brf", cache = TRUE, fig.cap = "Impact of plot length on profit estimated by BRF"----
g_plot_length_brf <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_brf_dif:=pi_brf-pi_opt,
         plot_length:=factor(plot_length,levels=c("12","24","48","Strip"))) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_brf_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=., aes(x=plot_length, y=mean, group=num_treatment)) +
  geom_line() +
  geom_point() +
  scale_x_discrete("plot_length", labels=c("12","24","48","Strip")) +
  facet_grid(num_treatment~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Plot Length') + 
  ylab('Estimated Deficit Profit \n by Balanced Random Forest (Kg/Ha)')
g_plot_length_brf


## ----fig.id = "plot-length-gam", cache = TRUE, fig.cap = "Impact of plot length on profit estimated by GAM"----
g_plot_length_gam <- all_sim_result %>% 
  mutate(plot_length:=ifelse(plot_length==288, 144, plot_length),
         pi_gam_dif:=pi_opt_gam-pi_opt_homo,
         plot_length:=factor(plot_length,levels=c("12","24","48","Strip"))) %>%
  .[yield_accuracy==0.9,] %>%
  .[,list(mean=mean(pi_gam_dif)), by=c('field_size', 'plot_length', 'num_treatment')]%>%
  ggplot(data=., aes(x=plot_length, y=mean, group=num_treatment)) +
  geom_line() +
  geom_point() +
  scale_x_discrete("plot_length", labels=c("12","24","48","Strip")) +
  facet_grid(num_treatment~field_size, labeller = as_labeller(var_names)) +
  theme(strip.text.y = element_text(angle = 360)) +
  scale_fill_discrete(name='Number of Treatments') + 
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Plot Length') + 
  ylab('Estimated Deficit Profit \n by Generalized Additive Model (Kg/Ha)')

g_plot_length_gam


## ----fig.id = "yield-accuracy-gwr", cache = TRUE, fig.cap = "Impact of yield monitor accuracy on profit estimated by GWR"----
g_yield_accuracy_gwr <- all_sim_result %>% 
  ggplot(data=.) +
  geom_density(aes(x=pi_gwr, fill=factor(yield_accuracy)),alpha=0.4) +
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Estimated Profit by Geographically Weighted Regression (Kg/Ha)') +
  ylab('Density') 

g_yield_accuracy_gwr


## ----fig.id = "yield-accuracy-brf", cache = FALSE, fig.cap = "Impact of yield monitor accuracy on profit estimated by BRF"----
g_yield_accuracy_brf <- all_sim_result %>% 
  ggplot(data=.) +
  geom_density(aes(x=pi_brf, fill=factor(yield_accuracy)),alpha=0.4) +
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Estimated Profit by Balanced Random Forest (Kg/Ha)') +
  ylab('Density') 

g_yield_accuracy_brf


## ----fig.id = "yield-accuracy-gam", cache = FALSE, fig.cap = "Impact of yield monitor accuracy on profit estimated by GAM"----
g_yield_accuracy_gam <- all_sim_result %>% 
  ggplot(data=.) +
  geom_density(aes(x=pi_opt_gam, fill=factor(yield_accuracy)),alpha=0.4) +
  theme(
    legend.position = 'bottom'
  ) +
  xlab('Estimated Profit by Generalized Additive Model (Kg/Ha)') +
  ylab('Density') 

g_yield_accuracy_gam

