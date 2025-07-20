############SPATIAL REGRESSION MODELS FOR AEREAL DATA###########################
#references
#David Leydet (2022-10-27) <https://rpubs.com/leydetd/spatialregressionI>
#David Leydet (2022-11-03) <https://rpubs.com/leydetd/spatialregressionII>
#Carlos Mendez (2020) <https://rpubs.com/quarcs-lab/tutorial-spatial-regression>
#Carlos Mendez (2020) <https://rpubs.com/quarcs-lab/spatial-autocorrelation>
#Data Science Book (2020)<https://github.com/CartoDB/data-science-book/blob/master/Chapter%201-2/Discrete%20Spatial%20Models.ipynb>
#<https://cran.r-project.org/web/packages/SDPDmod/vignettes/spatial_model.html>
#Maria R. Koldasheva and Nikolai A. Popov. (2023)<https://rpubs.com/Nick_Popov/thesis_spatial_dynamic>

# Load required libraries
library(sf) # Simple features for spatial data
library(spdep) # Spatial dependence tools
library(spatialreg) # Spatial regression models
library(ggplot2) # Plotting library
library(dplyr) # Data manipulation
library(gridExtra) # For arranging multiple plots
library(stats)
library(ggspatial) # For adding north arrow and scale to maps

# Set the path for saving plots
save_path <- "/home/edier/GoogleDrive/INVESTIGACION/PAPERS/ELABORACION/PAPER_SAR/Figures/"

# Load shapefile data (Area of Interest)
aoi = st_read("/home/edier/GoogleDrive/INVESTIGACION/PAPERS/ELABORACION/PAPER_SAR/Data/df_catchments_kmeans.gpkg", quiet = TRUE)

# Asegúrate de que la columna sea de tipo 'factor'
aoi$landcovermedian <- as.factor(aoi$landcovermedian)

# Change the reference level of 'landcovermedian' to a desired level, e.g., "grass"
aoi$landcovermedian <- relevel(aoi$landcovermedian, ref = "grass")

# Scale selected columns to standardize them for comparison
aoi2 <- aoi %>% mutate(across(c('rainfallAnnual_mean', 'elev_mean', 'rel_mean', 'area'), ~(scale(.) %>% as.vector)))

# Log transform the landslide variable to normalize the data
aoi2$logy = log(aoi$lands + 1)

# Store the geometry of the polygons for later use
aoi.geom = st_geometry(aoi)

# Calculate centroids of each polygon for neighbor identification
aoi.coords = st_centroid(aoi.geom)

# Create spatial neighbors using k-nearest neighbors (k = 5)
nb_k5 = knn2nb(knearneigh(aoi.coords, k = 5))

# Convert neighbors list to spatial weights
nb_k5_list = nb2listw(nb_k5)

# Run the Global Moran's I test to check for spatial autocorrelation
# Randomization assumption is TRUE as there are no known spatial trends
global_moran = moran.test(aoi2$logy, listw = nb_k5_list, alternative = "two.sided", randomisation = TRUE)
print(global_moran)

# Run Moran's I test multiple times using Monte Carlo simulations
# This helps determine if observed values are significantly different from random
moran_mc_simulation = moran.mc(aoi2$logy, listw = nb_k5_list, alternative = "greater", nsim = 999)
print(moran_mc_simulation)

# Create a Moran scatterplot to visualize the relationship between areas and their neighbors
moran.plot(aoi2$logy, nb_k5_list, labels = as.character(aoi$ID_CUENCA), xlab = "Landslide Density", ylab = "Lagged X")

# Fit a basic linear regression model
# The model predicts landslide log-density based on elevation, slope, and rainfall variables
col.fit1 = lm(logy ~ elev_mean + rel_mean + area, data = aoi2)
summary(col.fit1)

## morans I Monte Carlo simulations for RESIDUALS
moran.mc(residuals(col.fit1),nsim = 999,listw = nb_k5_list,alternative = "greater")

## The Lagrange multiplier test for RESIDUAL
lmt = lm.LMtests(col.fit1, nb_k5_list, test = c("LMerr", "LMlag"))
summary(lmt)

#Spatial Lag Model
col.fit2 = lagsarlm(logy ~ elev_mean + rel_mean + area, data = aoi2,listw = nb_k5_list)
summary(col.fit2,Nagelkerke=T)
#moran.mc(residuals(col.fit2),nsim = 999,listw = nb_k5_list,alternative = "greater")
#impacts(col.fit2, listw = nb_k5_list)
#summary(impacts(col.fit2, listw = nb_k5_list, R=500), zstats = TRUE) 

#Spatial Error Model
col.fit3 = errorsarlm(logy ~ elev_mean + rel_mean + area, data = aoi2,listw = nb_k5_list)
summary(col.fit3,Nagelkerke=T)
#impacts(col.fit3, listw = nb_k5_list)
#moran.mc(residuals(col.fit3),nsim = 999,listw = nb_k5_list,alternative = "greater")
#In the output of the function, note the value of lambda, the autoregressive coefficient representing the strength of autocorrelation in the residuals of a linear model.

#Spatial Durbin Lag Model
## Correlation between the dependent variable and the neighboring independent variables
## This also uses the lagsarlm() function, but with the parameter type set to ‘mixed’, to specify a Spatial Durbin lag model:
col.fit4 = lagsarlm(logy ~ elev_mean + rel_mean + area, data = aoi2,listw = nb_k5_list,type = "mixed")
summary(col.fit4,Nagelkerke=T)
#moran.mc(residuals(col.fit4),nsim = 999,listw = nb_k5_list,alternative = "greater")
#impacts(col.fit4, listw = nb_k5_list)
#summary(impacts(col.fit4, listw = nb_k5_list, R=500), zstats = TRUE) 

#SLX Spatially Lagged X
col.fit5 = lmSLX(logy ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = nb_k5_list)
summary(col.fit5,Nagelkerke=T)
#AIC(col.fit5)
#moran.mc(residuals(col.fit5),nsim = 999,listw = nb_k5_list,alternative = "greater")
#impacts(col.fit5, listw = nb_k5_list)
#summary(impacts(col.fit5, listw = nb_k5_list, R=500), zstats = TRUE) 

#Spatial Durbin Error
col.fit6 <- errorsarlm(logy ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = nb_k5_list, etype = "emixed")
summary(col.fit6,Nagelkerke=T)
#moran.mc(residuals(col.fit6),nsim = 999,listw = aoi.listw,alternative = "greater")

#Mansky (all inclusive - not recommended)
col.fit7 <- sacsarlm(logy ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = nb_k5_list, type="sacmixed") 
summary(col.fit7,Nagelkerke=T)
#moran.mc(residuals(col.fit7),nsim = 999,listw = nb_k5_list,alternative = "greater")
#impacts(col.fit7, listw = nb_k5_list)
#summary(impacts(col.fit7, listw = nb_k5_list, R=500), zstats = TRUE) 

#SARAR, Kelejian-Prucha, Cliff-Ord, or SAC
col.fit8 <- sacsarlm(logy ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = nb_k5_list, type="sac")
summary(col.fit8,Nagelkerke=T)
#moran.mc(residuals(col.fit8),nsim = 999,listw = nb_k5_list,alternative = "greater")
#summary(impacts(col.fit8, listw = nb_k5_list, R=500), zstats = TRUE) 

##################################

#regimes-Spatial Lag Model-cuenca
col.fit9 <- lagsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list)
summary(col.fit9,Nagelkerke=T)

#regimes-Spatial Lag Model-cluster
col.fit10 <- lagsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list)
summary(col.fit10,Nagelkerke=T)

#regimes-Spatial Error Model-cuenca
col.fit11 = errorsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list)
summary(col.fit11,Nagelkerke=T)

#regimes-Spatial Error Model-cluster
col.fit12 = errorsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list)
summary(col.fit12,Nagelkerke=T)

#regimes-Mansky-knn5
col.fit13 <- sacsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list, type="sacmixed") 
summary(col.fit13,Nagelkerke=T)

#regimes-Mansky-cuenca
col.fit14 <- sacsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list, type="sacmixed") 
summary(col.fit14,Nagelkerke=T)

#regimes-SAC-cuenca
col.fit15 <- sacsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list, type="sac")
summary(col.fit15)

#regimes-SAC-knn5
col.fit16 <- sacsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list, type="sac")
summary(col.fit16,adj.se=T,Nagelkerke=T)

#regimes-SLX-cuenca
col.fit17 = lmSLX(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list)
summary(col.fit17,Nagelkerke=T)
AIC(col.fit17)

#regimes-SLX-knn5
col.fit18 = lmSLX(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list)
summary(col.fit18,Nagelkerke=T)
AIC(col.fit18)

#regimes-SDM-cuenca
col.fit19 = lagsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = nb_k5_list,type = "mixed")
summary(col.fit19,Nagelkerke=T)

#regimes-SDM-knn5
col.fit19 = lagsarlm(logy ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = nb_k5_list,type = "mixed")
summary(col.fit19,Nagelkerke=T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##########FIGURE 5#####################
# Plot the geometry of the areas and their neighbor relationships and save the plot
neighbor_plot <- ggplot() +
  geom_sf(data = aoi.geom, fill = NA, color = "black") +
  geom_segment(data = do.call(rbind, lapply(1:length(nb_k5), function(i) {
    data.frame(
      x = st_coordinates(aoi.coords)[i, 1],
      y = st_coordinates(aoi.coords)[i, 2],
      xend = st_coordinates(aoi.coords)[nb_k5[[i]], 1],
      yend = st_coordinates(aoi.coords)[nb_k5[[i]], 2]
    )
  })), aes(x = x, y = y, xend = xend, yend = yend), color = "red", alpha = 0.5) +
  geom_sf(data = st_as_sf(aoi.coords), color = "blue", size = 1.5) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm")) +
  annotation_scale(location = "br", style = "ticks")

ggsave(filename = paste0(save_path, "../figFinal/Fig5.png"), plot = neighbor_plot, width = 10, height = 7, units = "in", dpi = 500)

############Fig 14###############

# --- 1. import data
if (!all(c("predicted_lands", "residuals") %in% names(aoi2))) {
  aoi2$predicted_logy <- as.vector(predict(col.fit3))
  aoi2$residuals <- as.vector(residuals(col.fit3))
  aoi2$predicted_lands <- exp(aoi2$predicted_logy) - 1
}

# --- 2. Create figures
g_predicted <- ggplot() +
  geom_sf(data = aoi2, aes(fill = predicted_lands), color = "black", linewidth = 0.1) +
  annotation_scale(location = "bl", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9 * 0.7, "cm"), width = unit(0.8 * 0.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  scale_fill_gradient(low = "white", high = "red", name = "Fitted landslide abundance") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.4), color = NA),
    
    # CAMBIO AQUÍ: Posición y estilo de la etiqueta (A, B, etc.)
    plot.tag.position = c(0.96, 0.05),
    plot.tag = element_text(face = "bold", size = 14)
  ) +
  # CAMBIO AQUÍ: Añadir la etiqueta "A"
  labs(x = NULL, y = NULL, tag = "A")


g_residuals <- ggplot() +
  geom_sf(data = aoi2, aes(fill = residuals), color = "black", linewidth = 0.1) +
  annotation_scale(location = "bl", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.9 * 0.7, "cm"), width = unit(0.8 * 0.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkgreen",
                       midpoint = 0, name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.4), color = NA),
    
    # CAMBIO AQUÍ: Posición y estilo de la etiqueta (A, B, etc.)
    plot.tag.position = c(0.96, 0.05),
    plot.tag = element_text(face = "bold", size = 16)
  ) +
  # CAMBIO AQUÍ: Añadir la etiqueta "B"
  labs(x = NULL, y = NULL, tag = "B")

# --- 4. Combine figures
final_plot_themed <- gridExtra::grid.arrange(g_predicted, g_residuals, nrow = 1)

ggsave(
  filename = paste0(save_path, "/figFinal/Fig14.png"),
  plot = final_plot_themed,
  width = 12,
  height = 6,
  units = "in",
  dpi = 500
)

