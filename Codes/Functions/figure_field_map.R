

#'/*=================================================*/
#'
#'/*                    Field Map				     */
#' 
#'/*=================================================*/

# raster map: plot the map when you first generate the data using gstat() function
cell_data <- gen_coefs(mean=200,psill=1000,range=sp_range,coef_name='N_star',nsim=B) %>%
    #>>> normalize <<<
    .[,sd_b:=sd(N_star),by=sim] %>%
    .[,mean_b:=mean(N_star),by=sim] %>%
    .[,p:=pnorm(N_star,mean=mean_b, sd=sd_b)] %>%
    .[,N_star:=100+p*150] %>%
    .[,c('cell_id','sim','X','Y','N_star')]
cell_data$N_star %>% hist(breaks=100)cell

rasterFromXYZ(cell_data[sim==1, c("X","Y","N_star")]) %>% plot()

{
  D <- cell_data[sim==1 & Y>500 & Y<505 & X>0 & X<7, c("X","Y","N_star")]
  # coordinates(D) = ~ X+Y
  R <- rasterFromXYZ(D)
}

# sf map: plot the map when you merge the cell level data to the field sf\

data <-
  cell_data[data.table(field), on = "cell_id"]

ggplot() +
    geom_sf(data = data, aes(fill = N_star), size=0.1, color=NA) +
    scale_fill_viridis_c(name="unit", direction = -1,
                         guide=guide_colorbar(frame.colour="black")) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position='right',
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=10),
        axis.text=element_text(color='black'))


#'/*===================================================*/
#'
#'/*                Field Layout				     */
#' 
#'/*=================================================*/

# after you get the field sf data ready
f<-field

#'/*--------------------------------------*/
#'#		Old-fashion plot()
#'/*--------------------------------------*/

#---sf to sp
f <- as_Spatial(f)
# save polygons
library(rgdal)
setwd(paste0(getwd(),"/Graph/designs/"))
writeOGR(f, dsn=".", layer="field_map", 
         driver="ESRI Shapefile", overwrite_layer=TRUE)
setwd(wd)
# -> BTW, writeOGR sucks in dsn setting; it only works for dsn="."

#---coloring plots as white & gray
f$even <- (f$plot_row_id + f$plot_col_id)%%2 == 0
f$color <- "white"
f$color[f$even] <- "grey90"
#---plot sp: manually draw plot polygon
{
    png(file='./Results/field_layout.png',
        width=7.2, height=4.32, units="in", res=300)
    par(mar=c(0,0,0,0))
    plot(f, border="grey60", col=f$color, lwd=0.1)
    
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



