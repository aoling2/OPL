#' title: "Graphs of trial design simulation results"
#'

rm(list=ls())

#/*=================================================*/
#' #                  Preparation
#/*=================================================*/

#=== Packages ===#
library(sf)
library(raster)
library(ggplot2)
library(data.table)
library(magrittr)
library(viridis)
library(dplyr)

#=== Set working directory ===#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
wd <- getwd()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#'/*================================================*/
#'
#'/*              Load and Prepare Data             */
#' 
#'/*================================================*/

#=== load results data ===#
all_results <- readRDS(
	file=paste0('./Results',
				'/base_1-1000_results_Sph_QP_pCpN_0.197_0.882.rds')) 


#====== data preparation ======#
pi_data <- all_results %>%
    .[order(range, fieldsize, design),] %>%
    #===keep columns
    dplyr::select(., pi_opt, pi_gwr, pi_rf_perfect,
                  design, range, fieldsize, sim) %>%
    #===relative profits
    .[, pi_gwr := pi_gwr-pi_opt] %>%
    .[, pi_rf_perfect := pi_rf_perfect-pi_opt] %>%
    .[, pi_opt:=NULL] %>%
    #=== Wide to Long: melt()
    melt(id.vars=c('range','fieldsize','design','sim')) %>%
    data.table() %>%
    #===generate label variables
    .[variable=='pi_gwr', application:='GWR'] %>%
    .[variable=='pi_rf_perfect', application:='Random Forest'] %>%
    .[fieldsize=='full', fieldha:=31.2] %>%
    .[fieldsize=='half', fieldha:=31.2/2] %>%
    .[fieldsize=='quarter', fieldha:=31.2/4] %>%
    .[,v_txt:=paste0('field = ',fieldha, ' (ha)')] %>%
    .[,v_txt:=factor(v_txt, levels=c('field = 31.2 (ha)',
                                     'field = 15.6 (ha)',
                                     'field = 7.8 (ha)'))] %>%
    .[,h_txt:=paste0('Range = ', range, ' (meter)')] %>%
    print()

#=== graphing labels ===#
N_apply <- c('GWR','Random Forest')
N_model <- c('pi_gwr','pi_rf_perfect')

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#'/*=================================================*/
#'
#' #              Mean Profits Bar Plot
#' 
#'/*=================================================*/

#=== average profit across simulations ===#
mean_data <- pi_data %>%
  #===take average across simulations
  .[, .(value=mean(value)), by=c("design","fieldsize","range","variable",
                                 "application","v_txt","h_txt")] %>%
	print()

value_ls <- -seq(0,-mean_data[,min(value)/5] %>% ceiling()*5,by=5)


#=== one scenario: range==600, field size=full
# rank designs by mean profit
design_level <- mean_data %>%
	.[range==600,] %>%
	.[fieldsize=='full',] %>%
	.[variable=='pi_gwr',] %>%
	.[order(-value),] %>%
	.$design

mean_data %>%
	.[range==600,] %>%
	.[fieldsize=='full',] %>%
	.[, design:=factor(design, levels=design_level)] %>%
  ggplot(data=., 
		   aes(x=application,y=value,fill=design)) +
	geom_bar(stat="identity", position=position_dodge2(reverse=TRUE)) +
  scale_x_discrete(position = "top") +
	coord_flip() +
	geom_text(aes(label=round(value,1)), 
			  position=position_dodge2(0.9,reverse=TRUE), hjust=-0.1, size=3)+
  geom_text(aes(x=application, y=0, label=design), 
            position=position_dodge2(0.9,reverse=TRUE), hjust=1.1, size=3) +
  ylab('Average Profit Relative to True Optimal ($/ha)') +
	xlab('') +
	scale_fill_discrete(name='Design') +
	scale_y_continuous(breaks = value_ls, label=value_ls) +
	theme(
		legend.position='none',
		legend.title = element_text(size=12),
		legend.key.size = unit(0.4, 'cm'),
		legend.text = element_text(size=10),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text=element_text(color='black')
	) 

ggsave(file=paste0('./Graph/mean profits/profits_combine_range600_full.png'),
	   height=6,width=6.5)


#=== all scenarios bar plot ===#
design_level_ls <- list(
  c("Latin Square","Alternate Block","Repeated Block","Checkerboard",
    "Randomized Block","Strip Fixed","Completely Random","Strip Random","Cascade Plot"),
  c("Alternate Block","Latin Square","Checkerboard","Repeated Block",
    "Randomized Block","Completely Random", "Strip Fixed","Strip Random","Cascade Plot")
)
for(i in 1:2){
  mean_data %>%
    .[application==N_apply[i],] %>%
    .[, design:=factor(design, levels=design_level_ls[[i]])] %>%
    ggplot(data=., 
         aes(x=application, y=value, fill=design)) +
    geom_bar(stat="identity", position=position_dodge2(reverse=TRUE)) +
    coord_flip() +
    facet_grid(v_txt~h_txt) +
    geom_text(aes(label=round(value,1)), 
              position=position_dodge2(0.9,reverse=TRUE), hjust=-0.1, size=3)+
    ggtitle(N_apply[i]) +
    ylab('Average Profit Relative to True Optimal VRA ($/ha)') +
    xlab('') +
    scale_fill_discrete(name='Design') +
    # scale_fill_discrete(guide=guide_legend(nrow=2,byrow=TRUE,title.vjust=0.5)) +
    theme(
      legend.position='right',
      legend.title = element_text(size=12),
      legend.key.size = unit(0.4, 'cm'),
      legend.text = element_text(size=10),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text=element_text(color='black')
    ) 
  ggsave(file=paste0('./Graph/mean profits/profits_',
                     N_apply[i],
                     '_combine','.png'),
         height=8,width=6.5)
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#/*=================================================*/
#'
#' #        Distribution of Simulated Profits
#' 
#/*=================================================*/

#=== one scenario: range==600, field size=full
pi_data_sc <- pi_data %>%
	.[range==600,] %>%
	.[fieldsize=='full',] %>%
  .[, linegroup:=ifelse(
    design%in%c("Alternate Block","Checkerboard","Repeated Block","Latin Square"),
    "High","Low")] %>%
  .[, design:=factor(design, levels=design_level)] %>%
  print()
for(i in 1:2){
	ggplot() +
		stat_density(data=pi_data_sc[variable==N_model[i],], 
					 aes(x=value, colour=design, linetype=linegroup),
					 position="identity", geom="line", size=0.5) +
		# geom_vline(data=mean_data, aes(xintercept=value, colour=design, linetype=design),
		# 		   size=1, show.legend=FALSE) +
		xlim(-75,-10) +
		xlab('Relative Profit to True Optimal ($/ha)') +
		ylab("Density") +
		labs(colour="Design", linetype="Group") +
		ggtitle(N_apply[i]) +
		theme(
		  legend.position='right',
		  legend.title = element_text(size=12),
		  legend.key.width = unit(1,"cm"),
		  legend.text = element_text(size=10),
		  plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		  axis.text=element_text(color='black')
		  )
		ggsave(file=paste0('./Graph/distribute/profits_kernel_',
					   N_apply[i],
					   '_range600_full.png'),
		   height=3,width=6.5)
}

# stacked figure
ggplot() +
  stat_density(data=pi_data_sc, 
               aes(x=value, colour=design, linetype=linegroup),
               position="identity", geom="line", size=0.5) +
  facet_wrap(~application, nrow=3) +
  xlim(-75,-10) +
  xlab('Relative Profit to True Optimal ($/ha)') +
  ylab("Density") +
  labs(colour="Design", linetype="Group") +
  theme(
    legend.position='right',
    legend.title = element_text(size=12),
    legend.key.width = unit(1,"cm"),
    legend.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text=element_text(color='black')
    )
ggsave(file=paste0('./Graph/distribute/profits_kernel_',
                   '_range600_full.png'),
       height=7.5,width=6.5)


#=== all scenarios:
for(i in 1:2){
  pi_data %>%
    .[application==N_apply[i],] %>%
    .[, design:=factor(design, levels=design_level_ls[[i]])] %>%
    ggplot(data=.) +
		stat_density(
					 aes(x=value, colour=design),
					 position="identity", geom="line", size=0.5) +
    facet_grid(v_txt~h_txt) +
    xlim(-100,-10) +
		xlab('Relative Profit to True Optimal ($/ha)') +
		ylab("Density") +
		labs(colour="Design", linetype="Design") +
		theme(
			legend.position='bottom',
			legend.key.width = unit(1,"cm"),
			legend.title = element_blank(),
			axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
			axis.text=element_text(color='black')
		) 
	ggsave(file=paste0('./Graph/distribute/profits_kernel_',
					   N_apply[i],
					   '_combine.png'),
		   height=7.5,width=6.5)
}


#=== extreme values ===#
for(i in 1:2){
	pi_data %>% 
		.[range==600,] %>%
    .[application==N_apply[i],] %>%
    .[, design:=factor(design, levels=design_level_ls[[i]])] %>%
    ggplot(data=., aes(x=value, colour=design, fill=design)) +
		geom_histogram(alpha=0.6, binwidth = 2.5) +
		facet_wrap(~design) +
		xlim(-100,-50) +
		xlab('Relative Profit to True Optimal ($/ha)') +
		ylab("Counts") +
		labs(colour="Design") +
		theme(
			legend.position='na',
			legend.title = element_blank(),
			axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
		) 
	ggsave(file=paste0('./Graph/distribute/profits_extreme_',
					   N_apply[i],
					   '_combine.png'),
		   height=3,width=6.5)
}

#=== percentage of extreme values ===#

#=== one scenario: range==600, field size=full
pi_data %>%
  .[range==600,] %>%
  .[fieldsize=='full',] %>%
  .[, .(count=sum(value<(-40)), nsim=length(sim)), 
    by=.(design, fieldsize, range, application)] %>%
  .[, extreme_percent:=round(count/nsim*100,1)] %>%
  #=== Long to Wide: dcast()
  dcast(design~application, value.var="extreme_percent") %>%
  data.table() %>%
  .[order(`GWR`),] %>%
  print() %>%
  saveRDS(.,
          paste0('./Graph/Tables/extreme_percent.rds'))



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#/*=================================================*/
#'
#' #              Pairwise Difference
#' 
#/*=================================================*/

#=== difference between pairs of designs ===#
diff_data <- pi_data %>%
	#=== Long to Wide: dcast()
	dcast(range+fieldsize+application+sim~design, value.var="value") %>%
	#=== pairwise difference
	.[, `Alternate Block - Latin Square` := `Alternate Block` - `Latin Square`] %>%
	.[, `Alternate Block - Checkerboard` := `Alternate Block` - `Checkerboard`] %>%
	.[, `Alternate Block - Randomized Block` := `Alternate Block` - `Randomized Block`] %>%
	.[, `Alternate Block - Completely Random` := `Alternate Block` - `Completely Random`] %>%
	.[, `Alternate Block - Strip Fixed` := `Alternate Block` - `Strip Fixed`] %>%
	.[, `Alternate Block - Strip Random` := `Alternate Block` - `Strip Random`] %>%
	.[, `Alternate Block - Cascade Plot` := `Alternate Block` - `Cascade Plot`] %>%
	.[, `Alternate Block`:=NULL] %>%
	.[, `Latin Square`:=NULL] %>%
	.[, `Checkerboard`:=NULL] %>%
	.[, `Randomized Block`:=NULL] %>%
	.[, `Cascade Plot`:=NULL] %>%
	.[, `Strip Fixed`:=NULL] %>%
	.[, `Strip Random`:=NULL] %>%
	.[, `Completely Random`:=NULL] %>%
  data.table() %>%
	#=== back to Long again
	melt(id.vars=c('range','fieldsize','application','sim')) %>%
	#=== scenario: range=600, field size=full
	.[range==600,] %>%
	.[fieldsize=='full',] %>%
	print()

#=== difference distribution
label_data <- data.table(
  application=N_apply,
  minx=c(-5, -20, -20),
  maxx=c(10, 40, 40),
  largerx=c(3, 3, 3),
  largery=c(0.15, 0.02, 0.02),
  lowerx=c(-4, -15, -15),
  lowery=c(0.15, 0.02, 0.02)
)
for(i in 1:2){
	#=== kernel curve data ===#
	den_data <- diff_data %>%
		.[application==N_apply[i],] %>%
		.[, .(x=density(value)$x,
			  y=density(value)$y),
		  by=.(range,fieldsize,application,variable)] %>%
		print()
	#=== percentage sum data ===#
	sum_data <- diff_data %>%
		.[application==N_apply[i],] %>%
		.[, .(larger=sum(value>0)/length(value)), 
		  by=.(application, fieldsize, range, variable)] %>%
		.[, larger_prct:=paste0(100*round(larger,2),'%')] %>%
		.[, lower_prct:=paste0(100*round(1-larger,2),'%')] %>%
		print()
	
	ggplot()+
		geom_line(data=den_data, aes(x,y)) +
		geom_area(data=subset(den_data,x>0), aes(x,y), fill="red", alpha=0.5) +
		geom_area(data=subset(den_data,x<0), aes(x,y), fill="blue", alpha=0.5) +
		facet_wrap(~variable, ncol=3) +
		geom_text(data=sum_data, aes(x=label_data[i,largerx], y=label_data[i,largery], 
		                             label=larger_prct), 
				  hjust=0, size=4) +
		geom_text(data=sum_data, aes(x=label_data[i,lowerx], y=label_data[i,lowery], 
		                             label=lower_prct), 
				  hjust=0, size=4) +
		xlim(label_data[i,minx],label_data[i,maxx]) +
		xlab('Profit Difference between Two Designs ($/ha)') +
		ylab('Density') +
		ggtitle(N_apply[i]) +
		theme(
		  plot.title = element_text(hjust = 0.5),
		  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		  axis.text=element_text(color='black')
		  )
	ggsave(file=paste0('./Graph/pairwise/profits_diff_',N_apply[i],'.png'),
		   height=7.5,width=7.5)
}
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#'/*=================================================*/
#'
#'/*        Percent of Best and Worst Designs        */
#' 
#'/*=================================================*/

#=== best/worse design in each simulation ===#
pi_data_2 <- pi_data  %>%
	#===find max profit
	.[, vmax:=max(value), by=.(range, fieldsize, application, sim)] %>%
	#===find min profit 
	.[, vmin:=min(value), by=.(range, fieldsize, application, sim)] %>%
	.[order(range, fieldsize, application, sim, design),] %>%
	print()
best_data <- pi_data_2 %>%
	.[value==vmax, .(range, fieldsize, application, sim, design)] %>%
	.[, perform:='Most Profitable'] %>%
	#===randomly select in ties
	.[sample(nrow(.)),] %>%
	unique(., by=c('range', 'fieldsize', 'application', 'sim')) %>%
	print()
worst_data <- pi_data_2 %>%
	.[value==vmin, .(range, fieldsize, application, sim, design)] %>%
	.[, perform:='Least Profitable'] %>%
	#===randomly select in ties
	.[sample(nrow(.)),] %>%
	unique(., by=c('range', 'fieldsize', 'application', 'sim')) %>%
	print()
best_worst_data <- rbind(best_data, worst_data) %>%
	.[order(range, fieldsize, application, sim),] %>%
	#===generate label variables
	.[, perform:=factor(perform, levels=c('Most Profitable','Least Profitable'))] %>%
	print()

#=== scenario: range=600, fieldsize=full
best_worst_data %>%
	.[range==600,] %>%
	.[fieldsize=='full',] %>%
  .[, design:=factor(design, levels=design_level)] %>%
  ggplot(data=., aes(x=application, fill=design)) +
	geom_bar(position = "fill", color='black', size=0.25) +
	# scale_fill_manual(values = c("white","grey60","grey30","grey0")) +
	facet_wrap(~perform) +
	scale_fill_discrete(name='Design') +
	xlab("") +
	ylab("Percent") +
	theme(
		legend.position='right',
		legend.key.size = unit(0.4, 'cm'),
		legend.text = element_text(size=10),
		axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text=element_text(color='black')
	) 
ggsave(file=paste0('./Graph/pairwise/percent_best_worst.png'),
	   height=7.5,width=6.5)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#'/*=================================================*/
#'
#'/*              Possible Explanations             */
#' 
#'/*=================================================*/

#=== All tentative explaining statistic measures ===#
exp_var_df <- data.table(
  variable=c('MSTmin', 'NB_gini', 'van_Es_var', 
             'moran_I', 'variation_N', 'GR'),
  varname=c('Evenness of distribution',
            'Neighbor unbalance',
            'Spatial unbalance',
            'Morans I',
            'Local variation',
            'Gradation'  ) ) %>%
  .[, varname:=factor(varname, levels=varname)]

#=== Explanation across designs ===#
source('./Codes/figure_explain_all_design.R')
figure_explain_all_design(raw_df=all_results, 
               exp_var_df=exp_var_df, 
               exp_figure_title="explain_all_designs_gwr")

#=== Explanation within design (Completely Random) ===#
source('./Codes/figure_explain_random.R')
random_df <- all_results[design=="Completely Random",]
figure_explain_random(raw_df=random_df, 
               exp_var_df=exp_var_df, 
               exp_figure_title="explain_random_gwr")

#=== Explanation within design (max 20 vs. min 20 Latin Square) ===#
source('./Codes/figure_explain_max_min.R')
Latin_results <- readRDS(
  file=paste0('./Results/GWR_1-1000_latin_rep_results_Sph_QP_pCpN_0.197_0.882.rds')) %>%
  setnames(., c("h_var_2d"), c("van_Es_var")) 
exp_var_df <- data.table(
  variable=c('MSTmin', 'van_Es_var', 
             'moran_I', 'GR'),
  varname=c('Evenness of distribution',
            'Spatial unbalance',
            'Morans I',
            'Gradation'  ) ) %>%
  .[, varname:=factor(varname, levels=varname)]
figure_explain_max_min(raw_df=Latin_results, 
                      exp_var_df=exp_var_df, 
                      exp_figure_title="explain_max_min_gwr")

{
  # check the local distribution of N trail rates
  # use a moving window method to check "local"
  
  #=== square window dimension
  celln <- 2*6
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  # accidental correlation
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #====== correlation N~coef ======#
  corr_N_coef_dt <- readRDS(
    file=paste0('./Results/corr_N_coef_dt_',celln,'.rds')) %>%
    .[design=='Strip', design:='Strip Fixed'] %>%
    .[design!="Repeated Block", ] %>%
    print()
  mean_data <- corr_N_coef_dt %>% 
    .[, .(corr_mean=mean(abs(corr_N_coef))), by=.(design)] %>%
    .[, design_label:=paste0(design,' (',round(corr_mean,3),')')] %>%
    .[order(corr_mean),] %>%
    .[, design_label:=factor(design_label, levels=design_label)] %>%
    print()
  corr_N_coef_dt <- corr_N_coef_dt %>%
    .[mean_data, on='design'] %>%
    .[, design:=factor(design, levels=mean_data$design)] %>%
    .[, linegroup:=ifelse(
      design%in%c("Alternate Block","Checkerboard","Latin Square"),
      "High","Low")] %>%
    print()
  ggplot() +
    stat_density(data=corr_N_coef_dt, 
                 aes(x=corr_N_coef, colour=design, linetype=linegroup),
                 position="identity", geom="line", size=0.5) +
    # labs(colour=expression(paste("Design  ( ", bar(paste("|", rho, "|")), " )")),
    # 	 linetype=expression(paste("Design  ( ", bar(paste("|", rho, "|")), " )"))) +
    # xlab("\u03c1(N Trial Rates, True Response Parameter)") +
    labs(colour="Design", linetype="Group") +
    xlab("Correlation between N Rates and True Response Parameter") +
    ylab("Density") +
    theme(
      legend.position='right',
      legend.title = element_text(size=12),
      legend.key.width = unit(1,"cm"),
      legend.text = element_text(size=10),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text=element_text(color='black')
    )
  ggsave(file=paste0('./Graph/explain/corr_N_coef_',celln,'cell.png'),
         height=4.25,width=6.5)
  
  #====== correlation N~error ======#	
  corr_N_error_dt <- readRDS(
    file=paste0('./Results/corr_N_error_dt_',celln,'.rds')) %>%
    .[design=='Strip', design:='Strip Fixed'] %>%
    .[design!="Repeated Block", ] %>%
    print()
  mean_data <- corr_N_error_dt %>% 
    .[, .(corr_mean=mean(abs(corr_N_error))), by=.(design)] %>%
    .[, design_label:=paste0(design,' (',round(corr_mean,3),')')] %>%
    .[order(corr_mean),] %>%
    .[, design_label:=factor(design_label, levels=design_label)] %>%
    print()
  corr_N_error_dt <- corr_N_error_dt %>%
    .[mean_data, on='design'] %>%
    .[, design:=factor(design, levels=mean_data$design)] %>%
    .[, linegroup:=ifelse(
      design%in%c("Alternate Block","Checkerboard","Latin Square"),
      "High","Low")] %>%
    print()
  ggplot() +
    stat_density(data=corr_N_error_dt, 
                 aes(x=corr_N_error, colour=design, linetype=linegroup),
                 position="identity", geom="line", size=0.5) +
    # labs(colour=expression(paste("Design  ( ", bar(paste("|", rho, "|")), " )")),
    # 	 linetype=expression(paste("Design  ( ", bar(paste("|", rho, "|")), " )"))) +
    # xlab("\u03c1(N Trial Rates, Yield Error Term)") +
    labs(colour="Design", linetype="Group") +
    xlab("Correlation between N Rates and True Yield Error") +
    ylab("Density") +
    theme(legend.key.width = unit(1.5,"cm"))
  ggsave(file=paste0('./Graph/explain/corr_N_error_',celln,'cell.png'),
         height=4.25,width=6.5)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  # local variation
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  variation_N_dt <- readRDS(
    file=paste0('./Results/variation_N_dt_',celln,'.rds')) %>%
    .[design=='Strip', design:='Strip Fixed'] %>%
    .[design!="Repeated Block", ] %>%
    print()
  mean_data <- variation_N_dt %>% 
    .[, .(sd_mean=mean(abs(variation_N))), by=.(design)] %>%
    .[, design_label:=paste0(design,' (',round(sd_mean,1),')')] %>%
    .[order(-sd_mean),] %>%
    .[, design_label:=factor(design_label, levels=design_label)] %>%
    print()
  variation_N_dt <- variation_N_dt %>%
    .[mean_data, on='design'] %>%
    .[, design:=factor(design, levels=mean_data$design)] %>%
    .[, linegroup:=ifelse(
      design%in%c("Alternate Block","Checkerboard","Latin Square"),
      "High","Low")] %>%
    print()
  ggplot() +
    stat_density(data=variation_N_dt, 
                 aes(x=variation_N, colour=design, linetype=linegroup),
                 position="identity", geom="line", size=0.5) +
    # labs(colour=expression(paste("Design  ( ", bar(sigma), " )")),
    # 	 linetype=expression(paste("Design  ( ", bar(sigma), " )"))) +
    # xlab("\u03c3(N Trial Rates)") +
    labs(colour="Design", linetype="Group") +
    xlab("Standard Deviation of N Rates") + 
    ylab("Density") +
    theme(legend.key.width = unit(1.5,"cm"))
  ggsave(file=paste0('./Graph/explain/variation_N_',celln,'cell.png'),
         height=4.25,width=6.5)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  # table
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  mean_variation_N_dt <- readRDS(
    file=paste0('./Results/variation_N_dt_',celln,'.rds')) %>%
    .[design=='Strip', design:='Strip Fixed'] %>%
    .[design!="Repeated Block", ] %>%
    .[, .(sd_N=mean(abs(variation_N))), by=.(design)] %>%
    .[, sd_N:=round(sd_N, 1)] %>%
    print()
  mean_corr_N_error_dt <- readRDS(
    file=paste0('./Results/corr_N_error_dt_',celln,'.rds')) %>%
    .[design=='Strip', design:='Strip Fixed'] %>%
    .[design!="Repeated Block", ] %>%
    .[, .(corr_N_error=mean(abs(corr_N_error))), by=.(design)] %>%
    .[, corr_N_error:=round(corr_N_error, 3)] %>%
    print()
  corr_tab <-  mean_variation_N_dt %>%
    .[mean_corr_N_error_dt, on='design'] %>%
    .[order(-sd_N),] %>%
    setnames(names(.), c('Design','sd(N)','corr(N, yield error)')) %>%
    print() %>%
    saveRDS(., paste0('./Graph/Tables/corr_tab.rds'))
}





#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#'/*=================================================*/
#'
#'/*              Field Layout Mapping				 */
#' 
#'/*=================================================*/

## load field layout data
all_fields <- readRDS('./Data/all_fields.rds')
field <- all_fields[["Randomized Block"]] %>%
	.[, c('cid','aunit_id','plot_id','block_id','buffer',
		  'plot_row_id','plot_col_id')]


#'/*--------------------------------------*/
#'#		Old-fashion plot()
#'/*--------------------------------------*/

#---sf to sp
field <- as_Spatial(field)
# save polygons
library(rgdal)
setwd(paste0(getwd(),"/Graph/designs/"))
writeOGR(field, dsn=".", layer="field_map", 
		 driver="ESRI Shapefile", overwrite_layer=TRUE)
setwd(wd)
# -> BTW, writeOGR sucks in dsn setting; it only works for dsn="."

#---coloring plots as white & gray
field$even <- (field$plot_row_id + field$plot_col_id)%%2 == 0
field$color <- "white"
field$color[field$even] <- "grey90"
#---plot sp: manually draw plot polygon
{
	png(file='./Graph/designs/field_layout.png',
		width=7.2, height=4.32, units="in", res=300)
	par(mar=c(0,0,0,0))
	plot(field, border="grey60", col=field$color, lwd=0.1)

	# cell label
	x_coord <- 72*4+18 + c(0, 6, 6, 0, 0)
	y_coord <- 420 + c(0, 0, 6, 6, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(72*4+27, 420, labels="cell", adj=c(0,0), cex=0.75)

	# subplot label
	x_coord <- 72*4+12 + c(0, 12, 12, 0, 0)
	y_coord <- 396 + c(0, 0, 18, 18, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(72*4+25, 396+6, labels="subplot", adj=c(0,0), cex=0.75)

	# plot label
	x_coord <- 72*5 + c(0, 72, 72, 0, 0)
	y_coord <- 396 + c(0, 0, 18, 18, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(72*5+24, 396+6, labels="plot", adj=c(0,0), cex=0.75)

	# block label
	x_coord <- 72*4 + c(0, 72*2, 72*2, 0, 0)
	y_coord <- 378 + c(0, 0, 54, 54, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=2, lty=1)
	text(72*6+12, 378+20, labels="block", adj=c(0,0), cex=1.25)

	dev.off()
}
#---plot sp
{
	png(file='./Graph/designs/field_layout.png',
		width=7.2, height=4.32, units="in", res=300)
	par(mar=c(0,0,0,0))
	
	# cell map
	plot(field, border="grey30", col=field$color, lwd=0.1)
	# cell label
	x_coord <- 78 + c(0, 6, 6, 0, 0)
	y_coord <- 420 + c(0, 0, 6, 6, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(78+9, 420-0, labels="cell", adj=c(0,0), cex=0.75)
	
	# subplot unit
	# sp.subplot <- aggregate(field[,c('cid','aunit_id')], by=list(field$aunit_id), 
	# 						FUN=mean) %>% .[, c('aunit_id')]
	sp.subplot[1,] %>% plot(., add=TRUE, border="red", lwd=1.5, lty=1)
	# # subplot label
	# x_coord <- 84 + c(0, 12, 12, 0, 0)
	# y_coord <- 396 + c(0, 0, 18, 18, 0)
	# p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	# sps = SpatialPolygons(list(p1))
	# plot(sps, add=T, border="black", lwd=1.5, lty=1)
	# text(84-6, 396-12, labels="subplot", adj=c(0,0), cex=0.75)
	# 
	# # plot label
	# x_coord <- 108 + c(0, 36, 36, 0, 0)
	# y_coord <- 396 + c(0, 0, 18, 18, 0)
	# p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	# sps = SpatialPolygons(list(p1))
	# plot(sps, add=T, border="black", lwd=1.5, lty=1)
	# text(108+6, 396+6, labels="plot", adj=c(0,0), cex=0.75)
	# 
	# # block label
	# x_coord <- 72 + c(0, 72, 72, 0, 0)
	# y_coord <- 378 + c(0, 0, 54, 54, 0)
	# p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	# sps = SpatialPolygons(list(p1))
	# plot(sps, add=T, border="black", lwd=2, lty=1)
	# text(72+12, 378-21, labels="block", adj=c(0,0), cex=1.25)
	
	dev.off()
}
#---plot cell map
{
	png(file='./Graph/designs/field_cell.png',
		width=7.2, height=4.32, units="in", res=300)
	par(mar=c(0,0,0,0))
	plot(field, border="grey30", lwd=0.1)
	
	# cell label
	x_coord <- 78 + c(0, 6, 6, 0, 0)
	y_coord <- 420 + c(0, 0, 6, 6, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(78+9, 420-0, labels="cell", adj=c(0,0), cex=1.5)
	
	dev.off()
} 

#---plot block map
{
	png(file='./Graph/designs/field_plot.png',
		width=7.2, height=4.32, units="in", res=300)
	par(mar=c(0,0,0,0))
	plot(field, border="grey30", lwd=0.1)
	
	# cell label
	x_coord <- 78 + c(0, 6, 6, 0, 0)
	y_coord <- 420 + c(0, 0, 6, 6, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(78+9, 420-0, labels="cell", adj=c(0,0), cex=0.75)
	
	# plot label
	x_coord <- 108 + c(0, 36, 36, 0, 0)
	y_coord <- 396 + c(0, 0, 18, 18, 0)
	p1 <- Polygons(list(Polygon(cbind(x_coord, y_coord))), 1)
	sps = SpatialPolygons(list(p1))
	plot(sps, add=T, border="black", lwd=1.5, lty=1)
	text(108+6, 396+6, labels="plot", adj=c(0,0), cex=0.75)
	
	dev.off()
} 


#'/*===========================================*/
all_fields <- readRDS('./Data/all_fields.rds')
au_sf <- readRDS('./Data/analysis_unit_sf.rds')

# check a specific design
field_dgn <- all_fields[["random_block"]]

#=== simple sf plot ===#
plot(field_dgn["block_id"])
plot(field_dgn["plot_in_block_id"])

#=== ggplot of ids ===#
ggplot(data=field_dgn) +
    geom_sf(aes()) +
    geom_sf_text(aes(label = block_id), size=3)

#=== overlay aunit and cells ===#
ggplot() +
    geom_sf(data=field_dgn, aes()) +
    geom_sf(data=au_sf, fill=NA, size=2) +
    geom_sf_text(data=field_dgn, aes(label = plot_id), size=2)

#=== the GWR-MC paper codes ===#
all_fields_MC <- readRDS('./Data/all_fields_MC.rds')
field_dgn_MC <- all_fields_MC[["6-cell"]]
au_sf_MC <- readRDS('./Data/analysis_unit_sf_MC.rds')
#=== overlay aunit and cells ===#
ggplot() +
  geom_sf(data=field_dgn_MC, aes()) +
  geom_sf(data=au_sf_MC, fill=NA, size=2) +
  geom_sf_text(data=field_dgn_MC, aes(label = plot_id), size=2)

#
field_dgn; field_dgn_MC
#
field <- all_fields[[temp_design]]
ggplot() +
  geom_sf(data=field, aes()) +
  geom_sf(data=au_sf, fill=NA, size=2) +
  geom_sf_text(data=field, aes(label = plot_id), size=2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# BE VERY CAREFUL OF R'S STUPID 'stringsAsFactors' SETTING !!!!!!!!!!
#
# if not explicitly specifying, the 'expand.grid' will enforce the
# conversion of characters into factors.

design_expand <- expand.grid(
  design="random_block") %>%
  data.table()
temp_design <- design_expand[1,design]
# here temp_design is forced into a factor. Those it still shows as
# random_block, it is treated as 1 (there is only one level in the design_expand)
# so the all_fields[[temp_design]] below actually returns all_fields[[1]] !!!
field <- all_fields[[temp_design]]
ggplot() +
  geom_sf(data=field, aes()) +
  geom_sf(data=au_sf, fill=NA, size=2) +
  geom_sf_text(data=field, aes(label = plot_in_block_id), size=2)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#' 
#' 
#' #/*=================================================*/
#' #' #             Check N trial layout
#' #/*=================================================*/
#' 
#' # see "N mapping" in the `check_single_sim_QP.R`
#' 
#' 
#' 
#' 
#' 
#' 
#' #
#' 
#' 
#' #/*----------------------------------*/
#' #' ## Parameter Spatial Distribution
#' #/*----------------------------------*/
#' 
#' temp_range <- 400
#' coef_data <- readRDS(file=paste0('./Data/simulated_coesf_QD_sprange_',temp_range,'.rds'))
#' field_sf <- readRDS('./Data/field.rds')
#' 
#' # take the i-th simulation
#' i=111
#' 
#' #=== b1 ===#
#' coef_data %>% 
#'     .[sim==i, c("X","Y","b1")] %>%
#'     rasterFromXYZ() %>%
#'     plot()
#' 
#' 
#' # left_join(field,data[,.(cid,opt_N)],by='cid') %>%
#' #   select(opt_N) %>%
#' #   plot()
#' 
#' #=== b2 ===#
#' coef_data %>% 
#'     .[sim==i, c("X","Y","b2")] %>%
#'     rasterFromXYZ() %>%
#'     plot()
#' 
#' #=== opt_N ===#
#' pCorn <- 0.197
#' pN <- 0.882
#' coef_data %>% 
#'     .[sim==i, ] %>%
#'     .[, opt_N:=(pN/pCorn-b1)/(2*b2)] %>%
#'     .[, c("X","Y","opt_N")] %>%
#'     rasterFromXYZ() %>%
#'     plot()
#' 
#' 
#' 
#' #/*=================================================*/
#' #' #             True yield response function
#' #/*=================================================*/

#=== quadratic-Plateau response
gen_yield_QP <- function(b0,b1,b2,Nk,N){
  yield <- (N<Nk)*(b0+b1*N+b2*N^2) + (N>=Nk)*(b0+b1*Nk+b2*Nk^2)
  return(yield)
}
N = 0:300
y = gen_yield_QP(b=5000, b1=120, b2=-0.3, Nk=180, N)/1000
png(file="./Graph/parameters/quadratic_plateau_curve.png",
        width=4.25, height=3.25, units="in", res=300)
par(mar=c(4,4,1,0)+0.1)
plot(y~N, type='l', lwd=2,
	 xlab='N rate (kg/ha)', ylab='corn yield (mg/ha)')
dev.off()


#' 
#' 
#' 
#' 
#' source('./Codes/functions_parameters_test.R')
#' 
#' #/*=================================================*/
#' #' # Preparation
#' #/*=================================================*/
#' 
#' #=== field ===#
#' all_fields <- readRDS('./Data/all_fields.rds')
#' 
#' #=== geo coordinates of cells ===#
#' xy <- st_centroid(all_fields[[1]]) %>%
#'   st_coordinates() %>%
#'   data.table() %>%
#'   .[,cid:=all_fields[[1]]$cid]
#' 
#' #=== analysis unit sf ===#
#' au_sf <- readRDS('./Data/analysis_unit_sf.rds')
#' 
#' xy <- st_coordinates(st_centroid(au_sf))
#' 
#' au_sf <- cbind(au_sf, xy)
#' 
#' #=== weights matrix ===#
#' X <- au_sf$X; Y <- au_sf$Y
#' n <- length(X)
#' dist <- matrix(NA,n,n)
#' for(k in 1:n){
#'   for(j in 1:n){
#'     dist[k,j] <- sqrt((X[k]-X[j])^2 + (Y[k]-Y[j])^2)
#'   }
#' }
#' weight_matrix <- 1/dist*as.numeric(dist<30)
#' diag(weight_matrix) <- 0
#' weight_matrix <- weight_matrix/rowSums(weight_matrix)
#' weight_list <- mat2listw(weight_matrix)
#' 
#' #=== subfield zones ===#
#' meanX <- (max(X)+min(X))/2
#' meanY <- (max(Y)+min(Y))/2
#' au_sf$zone <- NA
#' au_sf$zone[X<meanX&Y<meanY] <- 1
#' au_sf$zone[X>meanX&Y<meanY] <- 2
#' au_sf$zone[X<meanX&Y>meanY] <- 3
#' au_sf$zone[X>meanX&Y>meanY] <- 4
#' 
#' #=== load functions and parameters ===#
#' source('./Codes/functions_parameters_test.R')
#' source('./Codes/gen_field_ids.R')
#' 
#' #=== prices ===#
#' pCorn <- price_table[3,pCorn]
#' pN <- price_table[3,pN]
#' 
#' pCpN <- paste0(pCorn,'_',pN)
#' 
#' #/*=================================================*/
#' #' # Main Simulation
#' #/*=================================================*/
#' #=== number of iterations ===#
#' B <- 1
#' 
#' comp_list <- expand.grid(
#'   psill=psill_ls,
#'   design=field_config_table[,"design_name"],
#'   range=range_ls,
#'   stringsAsFactors = FALSE
#' ) %>%
#'   data.table()
#' 
#' i=1
#' comp_sim <- function(i, pCorn, pN) {
#'   
#'   temp_pars <- comp_list[i,]
#'   temp_psill <- temp_pars[,psill]
#'   temp_design <- temp_pars[,design]
#'   temp_range <- temp_pars[,range]
#'   
#'   #/*----------------------------------*/
#'   #' ## Load parameters
#'   #/*----------------------------------*/
#'   #=== coefficients ===#
#' 
#'   #=== error ===#
#'   m_error <- readRDS(paste0('./Data/m_error_sprange_',temp_range,'.rds')) %>%
#'     .[,c('cid','sim',paste0('m_error_',temp_psill)),with=FALSE] %>%
#'     setnames(names(.),c('cid','sim','m_error'))
#'   
#'   #=== define the field specification ===#
#'   field <- all_fields[[temp_design]]
#'   field_dt <- cbind(field,field %>% st_centroid %>% st_coordinates) %>%
#'     data.table()
  



#================================================================
# # check field N map
# data <- data.frame(all_fields[[1]]) %>% data.table()
# #...
# fsf <- all_fields[[1]]
# fsf <- merge(fsf, data)
# plot(fsf["N"])

#================================================================







