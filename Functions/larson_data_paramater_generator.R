# Prepare Field Parameters for Simulation

## Load functions and library
source(here("./Functions/Functions_for_Prep.R"))
library("rgdal")
library("sp")
# library("raster")
# library("tmap")
# library("maptools")
library("dplyr")
# library("Hmisc")
library("sf")
# library("SpatialTools")
# library("geoR")
# library("scatterplot3d")
# library("spdep")
# library("rdist")
# library("fields")
# library("IPEC")
library("here")

## Read Cleaned Data###
field <- read_sf(here("Data","larson_eb2_data_cleaned_201760_10.shp"))%>%
	mutate(yield:=mean_y*1000, # convert yield to kg/ha
		n:=(mean_s*11.06*0.32+192)*1.12085, # convert nitrogen to kg/ha
		id:=1:nrow(.)) %>%
	select(id, yield, n)

lm_reg <- lm(yield~n+I(n^2), data = field)
summary(field$yield)
b_0 <- lm_reg$coefficients[1]
b_1 <- lm_reg$coefficients[2]
b_2 <- lm_reg$coefficients[3]
(pN/pCorn-b_1)/(2*b_2)
N_eval_ls<-seq(min(ceiling(field$n)),max(ceiling(field$n)), by = 0.1)%>%
	data.frame()


ggplot(N_eval_ls,aes(seq(min(ceiling(field$n)),max(ceiling(field$n)), by = 0.1)))+
	stat_function(fun=function(x) b_0+b_1*x+b_2*x^2)

N_eval_ls<-seq(1, 500, by = 0.1)%>%
	data.frame()

ggplot(N_eval_ls,aes(seq(1,500, by = 0.1)))+
	stat_function(fun=function(x) b_0+b_1*x+b_2*x^2)

Nk <- -b_1/(2*b_2)
Ymax <- b_0+b_1*Nk+b_2*Nk^2

saveRDS(lm_reg$coefficients, file = "lm_res")