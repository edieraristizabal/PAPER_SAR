############SPATIAL REGRESSION MODELS FOR AEREAL DATA###########################
#references
#David Leydet (2022-10-27) <https://rpubs.com/leydetd/spatialregressionI>
#David Leydet (2022-11-03) <https://rpubs.com/leydetd/spatialregressionII>
#Carlos Mendez (2020) <https://rpubs.com/quarcs-lab/tutorial-spatial-regression>
#Carlos Mendez (2020) <https://rpubs.com/quarcs-lab/spatial-autocorrelation>
#Data Science Book (2020)<https://github.com/CartoDB/data-science-book/blob/master/Chapter%201-2/Discrete%20Spatial%20Models.ipynb>
#<https://cran.r-project.org/web/packages/SDPDmod/vignettes/spatial_model.html>
#Maria R. Koldasheva and Nikolai A. Popov. (2023)<https://rpubs.com/Nick_Popov/thesis_spatial_dynamic>

library(RColorBrewer) #color palettes
library(sf) #simple features for spatial data
library(spdep)
library(spatialreg)
library(terra) # to support tmap
library(tmap) #mapping package
library(kableExtra) #table modification
library(ggplot2)
library(dplyr)
library(ggspatial) #para ponerle escala

aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi2 <- aoi %>% mutate_at(c('rainfallAnnual_mean','elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi2$lands_dens=aoi$lands_rec/aoi$Área

## Store the geometry of the polygons
aoi.geom = st_geometry(aoi)

## Store the centroids of the polygons
aoi.coords = st_centroid(aoi.geom)

## Here each region is linked to its k nearest neighbors irrespective of the distance between them
sy5_nb = knn2nb(knearneigh(aoi.coords, k = 5))
aoi.listw = nb2listw(sy5_nb)

plot(aoi.geom,reset = FALSE)
plot(sy5_nb, aoi.coords,add = TRUE,lwd=.2, col="blue", cex = .5)

## Run the Global Moran's I test to look for auto correlation. Randomization assumtion is set to TRUE as we don't know about any larger spatial trends in the data. Two sided meaning postive or negative autocorrelation. 
moran.test(aoi2$lands_dens,listw = aoi.listw,alternative = "two.sided",randomisation = TRUE)

## Run the test multiple times to compare our observed statistic compared to a random sampling of values across the dataset. 
## Once done, we look at the rank of the observed version of Moran’s I against those obtained from random resampling. If we are either at the high or low end of all these random realizations, it is highly likely that the observed distribution is significantly autocorrelated. As we are using a rank, we cannot use a two-sided test, but must specify if we believe the autocorrelation to be positive (‘greater’) or negative (‘less’). The number of simulations to run is given by the parameter nsim. Increasing this, will increase the precision of the p-value obtained (but take longer to run):
moran.mc(aoi2$lands_dens,listw = aoi.listw,alternative = "greater",nsim = 999)

## Now make the Moran scatterplot. This shows the relationship between the value of any given area and its neighbors. The slope of the fitted line is the value of Moran’s I:
moran.plot(aoi2$lands_dens,aoi.listw,labels = as.character(aoi$ID_CUENCA),xlab = "Landslide density",ylab = "Lagged X")

## We will now calculate the local Moran’s I to examine the spatial pattern of autocorrelation. This returns a statistic (and z-score) for each area considered.
lm1 = localmoran(aoi2$lands_dens,listw = aoi.listw,alternative = "two.sided")
head(lm1)

## Basic linear model
col.fit1 = lm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2)
summary(col.fit1)

## Use morans I Monte Carlo simulations
moran.mc(residuals(col.fit1),nsim = 999,listw = aoi.listw,alternative = "greater")

## The Lagrange multiplier test is used to assess whether the autocorrelation is in the values of the dependent variable or in its errors, and helps in the choice of which spatial regression model to use. 
## LMerr is testing the errors, LMlag is testing the dependent variable and a "lag" effect (spillover)
lmt = lm.LMtests(col.fit1, aoi.listw, test = c("LMerr", "LMlag"))
summary(lmt)

## Use the robust test to decide which of two is the more likely source for autocorrelation
#lmt_robust = lm.RStests(col.fit1, aoi.listw, test="all")
#summary(lmt_robust)

#Spatial Lag Model
col.fit2 = lagsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw)
summary(col.fit2,Nagelkerke=T)
#moran.mc(residuals(col.fit2),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit2, listw = aoi.listw)
#summary(impacts(col.fit2, listw = aoi.listw, R=500), zstats = TRUE) 

#Spatial Error Model
col.fit3 = errorsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw)
summary(col.fit3,Nagelkerke=T)
#moran.mc(residuals(col.fit3),nsim = 999,listw = aoi.listw,alternative = "greater")
#In the output of the function, note the value of lambda, the autoregressive coefficient representing the strength of autocorrelation in the residuals of a linear model.

#Spatial Durbin Lag Model
## Correlation between the dependent variable and the neighboring independent variables
## This also uses the lagsarlm() function, but with the parameter type set to ‘mixed’, to specify a Spatial Durbin lag model:
col.fit4 = lagsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw,type = "mixed")
summary(col.fit4,Nagelkerke=T)
#moran.mc(residuals(col.fit4),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit4, listw = aoi.listw)
#summary(impacts(col.fit4, listw = aoi.listw, R=500), zstats = TRUE) 


#SLX Spatially Lagged X
col.fit5 = lmSLX(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw)
summary(col.fit5,Nagelkerke=T)
#moran.mc(residuals(col.fit5),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit5, listw = aoi.listw)
#summary(impacts(col.fit5, listw = aoi.listw, R=500), zstats = TRUE) 

#Spatial Durbin Error
col.fit6 <- errorsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, etype = "emixed")
summary(col.fit6,Nagelkerke=T)
#moran.mc(residuals(col.fit6),nsim = 999,listw = aoi.listw,alternative = "greater")

#Mansky (all inclusive - not recommended)
col.fit7 <- sacsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit7,Nagelkerke=T)
moran.mc(residuals(col.fit7),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit7, listw = aoi.listw)
summary(impacts(col.fit7, listw = aoi.listw, R=500), zstats = TRUE) 

#SARAR, Kelejian-Prucha, Cliff-Ord, or SAC
col.fit8 <- sacsarlm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit8,Nagelkerke=T)
moran.mc(residuals(col.fit8),nsim = 999,listw = aoi.listw,alternative = "greater")
summary(impacts(col.fit8, listw = aoi.listw, R=500), zstats = TRUE) 


##################################

#regimes-Spatial Lag Model-NOMZH
col.fit9 <- lagsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(NOMZH), data = aoi2,listw = aoi.listw)
summary(col.fit9,Nagelkerke=T)

#regimes-Spatial Lag Model-cluster
col.fit10 <- lagsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit10,Nagelkerke=T)

#regimes-Spatial Error Model-NOMZH
col.fit11 = errorsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(NOMZH), data = aoi2,listw = aoi.listw)
summary(col.fit11,Nagelkerke=T)

#regimes-Spatial Error Model-cluster
col.fit12 = errorsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit12,Nagelkerke=T)

#regimes-Mansky-knn5
col.fit13 <- sacsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(knn5), data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit13,Nagelkerke=T)

#regimes-Mansky-NOMZH
col.fit14 <- sacsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(NOMZH), data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit14,Nagelkerke=T)

#regimes-SAC-NOMZH
col.fit15 <- sacsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(NOMZH), data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit15)

#regimes-SAC-knn5
col.fit16 <- sacsarlm(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(knn5), data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit16,adj.se=T,Nagelkerke=T)

#regimes-SLX-NOMZH
col.fit17 = lmSLX(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(NOMZH), data = aoi2,listw = aoi.listw)
summary(col.fit17,Nagelkerke=T)
AIC(col.fit17)

#regimes-SLX-knn5
col.fit18 = lmSLX(lands_dens ~ 0 + (elev_mean + slope_mean + RainfallDaysmean):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit18,Nagelkerke=T)
AIC(col.fit18)

#Comparing
aic.tbl = AIC(col.fit2, col.fit3, col.fit4,col.fit5,col.fit6,col.fit7,col.fit8)
rownames(aic.tbl) = c("Lag Model", "Error Model", "Durbin Model")
aic.tbl %>%kbl(caption = "AIC Comparison") %>%kable_classic_2(full_width = F, html_font = "arial")

aoi$res7=residuals(col.fit7)
aoi$fit7=fitted(col.fit7)

ggplot() + geom_sf(data=aoi,aes(fill=res7),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm"),) +
  scale_fill_gradientn(colors=c("white","orange"),name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=6),
    legend.title.align = 0,
    legend.position = c(0.3,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=10, vjust = .8, hjust = .5))


col.fit19=spautolm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, family="CAR")
summary(col.fit19,Nagelkerke=T)

col.fit20=spautolm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, family="SAR")
summary(col.fit20,Nagelkerke=T)
print(col.fit19)
