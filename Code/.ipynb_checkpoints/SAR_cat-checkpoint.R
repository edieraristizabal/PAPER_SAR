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
save_path <- "G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/PAPER_SAR/Figures/"

# Load shapefile data (Area of Interest)
aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg", quiet = TRUE)

# Change the reference level of 'landcovermedian' to a desired level, e.g., "grass"
aoi$landcovermedian <- relevel(aoi$landcovermedian, ref = "grass")

# Scale selected columns to standardize them for comparison
aoi2 <- aoi %>% mutate(across(c('rainfallAnnual_mean', 'elev_mean', 'slope_mean', 'RainfallDaysmean'), ~(scale(.) %>% as.vector)))

# Add a new column called 'y_density' with the rate of lands_rec/area
aoi2$y_density = aoi$lands_rec / aoi$area

# Log transform the landslide variable to normalize the data
aoi2$logy = log(aoi$lands + 1)

# Rank-based inverse normal transformation function
rank_inverse_normal <- function(x) {
  ranks <- rank(x, ties.method = "average")
  n <- length(ranks)
  qnorm((ranks - 0.5) / n)
}

# Apply the transformation to your 'y_density' column
aoi2$y_transf <- rank_inverse_normal(aoi2$y_density)

# Shift the transformed values to ensure all values are positive
shift_value <- -min(aoi2$y_transf) + 1 
aoi2$y_tp <- aoi2$y_transf + shift_value  # Apply the shift


# Plot histograms to visualize the distribution of the new variables
max_y <- max(table(cut(aoi2$y_density, breaks = 30)), table(cut(aoi2$logy, breaks = 30)), table(cut(aoi2$y_tp, breaks = 30)))

# Plot histograms to visualize the distribution of the new variables
p1 <- ggplot(aoi2, aes(x = y_transf)) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) +
  xlab("Density") +
  ylim(0, 250) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p2 <- ggplot(aoi2, aes(x = y_transf)) +
  geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
  xlab("Transformed Landslide Density") +
  ylim(0, 250) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p3 <- ggplot(aoi2, aes(x = y_tp)) +
  geom_histogram(fill = "lightcoral", color = "black", bins = 30) +
  xlab("Transformed Density") +
  ylim(0, 250) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Arrange all three plots in a 1-row, 3-column layout and save the output
combined_plot <- grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
ggsave(filename = paste0(save_path, "Histograms.png"), plot = combined_plot, width = 15, height = 5, units = "in", dpi = 500)

# Store the geometry of the polygons for later use
aoi.geom = st_geometry(aoi)

# Calculate centroids of each polygon for neighbor identification
aoi.coords = st_centroid(aoi.geom)

# Create spatial neighbors using k-nearest neighbors (k = 5)
nb_k5 = knn2nb(knearneigh(aoi.coords, k = 5))

# Convert neighbors list to spatial weights
nb_k5_list = nb2listw(nb_k5)

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

ggsave(filename = paste0(save_path, "AOI_Geometry_Neighbor_Connections.png"), plot = neighbor_plot, width = 10, height = 7, units = "in", dpi = 500)

# Run the Global Moran's I test to check for spatial autocorrelation
# Randomization assumption is TRUE as there are no known spatial trends
global_moran = moran.test(aoi2$y_density, listw = nb_k5_list, alternative = "two.sided", randomisation = TRUE)
print(global_moran)

# Run Moran's I test multiple times using Monte Carlo simulations
# This helps determine if observed values are significantly different from random
moran_mc_simulation = moran.mc(aoi2$y_density, listw = nb_k5_list, alternative = "greater", nsim = 999)
print(moran_mc_simulation)

# Create a Moran scatterplot to visualize the relationship between areas and their neighbors
moran.plot(aoi2$y_density, nb_k5_list, labels = as.character(aoi$ID_CUENCA), xlab = "Landslide Density", ylab = "Lagged X")

# Calculate local Moran's I to examine local spatial autocorrelation
# This provides a statistic for each individual area to identify clusters
local_moran = localmoran(aoi2$y_density, listw = nb_k5_list, alternative = "two.sided")
head(local_moran)

# Fit a basic linear regression model
# The model predicts landslide log-density based on elevation, slope, and rainfall variables
col.fit1 = lm(y_tp ~ elev_mean + slope_mean + landcovermedian, data = aoi2)
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
col.fit2 = lagsarlm(y_tp ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw)
summary(col.fit2,Nagelkerke=T)
#moran.mc(residuals(col.fit2),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit2, listw = aoi.listw)
#summary(impacts(col.fit2, listw = aoi.listw, R=500), zstats = TRUE) 

#Spatial Error Model
col.fit3 = errorsarlm(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw)
summary(col.fit3,Nagelkerke=T)
#moran.mc(residuals(col.fit3),nsim = 999,listw = aoi.listw,alternative = "greater")
#In the output of the function, note the value of lambda, the autoregressive coefficient representing the strength of autocorrelation in the residuals of a linear model.

#Spatial Durbin Lag Model
## Correlation between the dependent variable and the neighboring independent variables
## This also uses the lagsarlm() function, but with the parameter type set to ‘mixed’, to specify a Spatial Durbin lag model:
col.fit4 = lagsarlm(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw,type = "mixed")
summary(col.fit4,Nagelkerke=T)
#moran.mc(residuals(col.fit4),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit4, listw = aoi.listw)
#summary(impacts(col.fit4, listw = aoi.listw, R=500), zstats = TRUE) 

#SLX Spatially Lagged X
col.fit5 = lmSLX(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw)
summary(col.fit5,Nagelkerke=T)
AIC(col.fit5)
#moran.mc(residuals(col.fit5),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit5, listw = aoi.listw)
#summary(impacts(col.fit5, listw = aoi.listw, R=500), zstats = TRUE) 

#Spatial Durbin Error
col.fit6 <- errorsarlm(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw, etype = "emixed")
summary(col.fit6,Nagelkerke=T)
#moran.mc(residuals(col.fit6),nsim = 999,listw = aoi.listw,alternative = "greater")

#Mansky (all inclusive - not recommended)
col.fit7 <- sacsarlm(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit7,Nagelkerke=T)
moran.mc(residuals(col.fit7),nsim = 999,listw = aoi.listw,alternative = "greater")
#impacts(col.fit7, listw = aoi.listw)
summary(impacts(col.fit7, listw = aoi.listw, R=500), zstats = TRUE) 

#SARAR, Kelejian-Prucha, Cliff-Ord, or SAC
col.fit8 <- sacsarlm(y_transf ~ elev_mean + slope_mean + landcovermedian, data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit8,Nagelkerke=T)
moran.mc(residuals(col.fit8),nsim = 999,listw = aoi.listw,alternative = "greater")
summary(impacts(col.fit8, listw = aoi.listw, R=500), zstats = TRUE) 


##################################

#regimes-Spatial Lag Model-cuenca
col.fit9 <- lagsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw)
summary(col.fit9,Nagelkerke=T)

#regimes-Spatial Lag Model-cluster
col.fit10 <- lagsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit10,Nagelkerke=T)

#regimes-Spatial Error Model-cuenca
col.fit11 = errorsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw)
summary(col.fit11,Nagelkerke=T)

#regimes-Spatial Error Model-cluster
col.fit12 = errorsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit12,Nagelkerke=T)

#regimes-Mansky-knn5
col.fit13 <- sacsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit13,Nagelkerke=T)

#regimes-Mansky-cuenca
col.fit14 <- sacsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw, type="sacmixed") 
summary(col.fit14,Nagelkerke=T)

#regimes-SAC-cuenca
col.fit15 <- sacsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit15)

#regimes-SAC-knn5
col.fit16 <- sacsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw, type="sac")
summary(col.fit16,adj.se=T,Nagelkerke=T)

#regimes-SLX-cuenca
col.fit17 = lmSLX(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw)
summary(col.fit17,Nagelkerke=T)
AIC(col.fit17)

#regimes-SLX-knn5
col.fit18 = lmSLX(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw)
summary(col.fit18,Nagelkerke=T)
AIC(col.fit18)

#regimes-SDM-cuenca
col.fit19 = lagsarlm(y_transf ~ 0 + (elev_mean + slope_mean + landcovermedian):(cuenca), data = aoi2,listw = aoi.listw,type = "mixed")
summary(col.fit19,Nagelkerke=T)

#regimes-SDM-knn5
col.fit19 = lagsarlm(y_tp ~ 0 + (elev_mean + slope_mean + landcovermedian):(knn5), data = aoi2,listw = aoi.listw,type = "mixed")
summary(col.fit19,Nagelkerke=T)

#Comparing
aic.tbl = AIC(col.fit2, col.fit3, col.fit4,col.fit5,col.fit6,col.fit7,col.fit8)
rownames(aic.tbl) = c("Lag Model", "Error Model", "Durbin Model")
aic.tbl %>%kbl(caption = "AIC Comparison") %>%kable_classic_2(full_width = F, html_font = "arial")

aoi$res19=residuals(col.fit19)

# Extract the fitted values from the model in the transformed space
aoi$fit19_transformed <- fitted(col.fit19)

# Apply inverse of the rank inverse normal transformation (approximate to original scale)
# Since it's complex to exactly invert the rank-based normal transformation,
# you can map the predicted ranks back to the original values:
original_order <- order(rank(aoi$lands_rec / aoi$area, ties.method = "average"))
aoi$fit19_approx <- sort(aoi$fit19_transformed)[original_order]

# Plot residuals and approximated fitted values in a 1-row, 2-column layout
g1 <- ggplot() +
  geom_sf(data = aoi, aes(fill = res19), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_viridis_c(name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

g2 <- ggplot() +
  geom_sf(data = aoi, aes(fill = fit19_approx), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Fitted Values (Approx)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

combined_plot_2 <- grid.arrange(g1, g2, nrow = 1, ncol = 2)
ggsave(filename = paste0(save_path, "Residuals_and_Fitted_Values.png"), plot = combined_plot_2, width = 15, height = 7, units = "in", dpi = 500)

# Plot y_density and approximated fitted values in a 1-row, 2-column layout
g1 <- ggplot() +
  geom_sf(data = aoi2, aes(fill = y_density), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Y Density") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

g2 <- ggplot() +
  geom_sf(data = aoi, aes(fill = fit19_transformed), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Fitted Values (Approx)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

combined_plot_2 <- grid.arrange(g1, g2, nrow = 1, ncol = 2)
ggsave(filename = paste0(save_path, "Y_Density_and_Fitted_transf_Values.png"), plot = combined_plot_2, width = 15, height = 7, units = "in", dpi = 500)


########################################################################



col.fit19=spautolm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, family="CAR")
summary(col.fit19,Nagelkerke=T)

col.fit20=spautolm(lands_dens ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2,listw = aoi.listw, family="SAR")
summary(col.fit20,Nagelkerke=T)
print(col.fit19)
