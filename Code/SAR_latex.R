# Load required libraries
library(spatialreg)
library(spdep)
library(sf)
library(dplyr)
library(xtable)

# Load shapefile data (Area of Interest)
aoi = st_read("/home/edier/GoogleDrive/INVESTIGACION/PAPERS/ELABORACION/PAPER_SAR/Data/df_catchments_kmeans.gpkg", quiet = TRUE)

# Ensure the continuous variables (area, elev_mean, rel_mean) are numeric before scaling
aoi <- aoi %>%
  mutate(across(c(area, elev_mean, rel_mean), as.numeric))

# Standardize the continuous variables
aoi <- aoi %>%
  mutate(across(c(area, elev_mean, rel_mean), scale))

# Log transform the landslide variable to normalize the data
aoi$y_log = log(aoi$lands_rec + 1)

# Store the geometry of the polygons for later use
aoi.geom = st_geometry(aoi)

# Calculate centroids of each polygon for neighbor identification
aoi.coords = st_centroid(aoi.geom)

# Create spatial neighbors using k-nearest neighbors (k = 5)
nb_k5 = knn2nb(knearneigh(aoi.coords, k = 5))

# Convert neighbors list to spatial weights
nb_k5_list = nb2listw(nb_k5)

# Initialize an empty data frame to store the results
results <- data.frame(
  Parameter = c("Constant", "A_Coefficient", "E_Coefficient", "H_Coefficient", "Wx_A", "Wx_E", "Wx_H", "Rho", "Lambda", "Nagelkerke_R2", "AIC"),
  SAR = NA, SEM = NA, SAC = NA, SLX = NA, SDEM = NA, SDM = NA, GNS = NA,
  stringsAsFactors = FALSE
)

# Function to extract model summary information with standard errors and p-values
extract_model_summary <- function(model) {
  summary_model <- summary(model, Nagelkerke = TRUE)
  
  coefficients <- coef(summary_model)
  std_errors <- summary_model$Coef[, "Std. Error"]  # Extract standard errors
  p_values <- summary_model$Coef[, "Pr(>|z|)"]  # Extract p-values
  rho <- if (!is.null(summary_model$rho)) summary_model$rho else NA
  lambda <- if (!is.null(summary_model$lambda)) summary_model$lambda else NA
  nagelkerke <- if (!is.null(summary_model$NK)) summary_model$NK else NA
  aic <- AIC(model)
  
  return(list(coefficients = coefficients, std_errors = std_errors, p_values = p_values, rho = rho, lambda = lambda, nagelkerke = nagelkerke, aic = aic))
}

# Format coefficients with standard errors and bold significant ones using sprintf for precise formatting
format_coef <- function(coef, std_err, p_value) {
  if (!is.na(p_value) && p_value < 0.05) {
    return(paste0("\\textbf{", sprintf("%.3f", coef), "} (", sprintf("%.2f", std_err), ")"))
  } else {
    return(paste0(sprintf("%.3f", coef), " (", sprintf("%.2f", std_err), ")"))
  }
}

# Spatial Lag Model (SAR)
SAR_model <- lagsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list)
SAR_summary <- extract_model_summary(SAR_model)

results$SAR <- c(
  format_coef(SAR_summary$coefficients[1], SAR_summary$std_errors[1], SAR_summary$p_values[1]),  # Constant
  format_coef(SAR_summary$coefficients[3], SAR_summary$std_errors[3], SAR_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(SAR_summary$coefficients[2], SAR_summary$std_errors[2], SAR_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(SAR_summary$coefficients[4], SAR_summary$std_errors[4], SAR_summary$p_values[4]),  # H_Coefficient (rel_mean)
  NA, NA, NA,  # Wx_A, Wx_E, Wx_H (Not applicable for SAR)
  round(SAR_summary$rho, 3),  # Rho
  NA,  # Lambda (Not applicable for SAR)
  round(SAR_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SAR_summary$aic, 3)  # AIC
)

# Spatial Error Model (SEM)
SEM_model <- errorsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list)
SEM_summary <- extract_model_summary(SEM_model)

results$SEM <- c(
  format_coef(SEM_summary$coefficients[1], SEM_summary$std_errors[1], SEM_summary$p_values[1]),  # Constant
  format_coef(SEM_summary$coefficients[3], SEM_summary$std_errors[3], SEM_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(SEM_summary$coefficients[2], SEM_summary$std_errors[2], SEM_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(SEM_summary$coefficients[4], SEM_summary$std_errors[4], SEM_summary$p_values[4]),  # H_Coefficient (rel_mean)
  NA, NA, NA,  # Wx_A, Wx_E, Wx_H (Not applicable for SEM)
  NA,  # Rho (Not applicable for SEM)
  round(SEM_summary$lambda, 3),  # Lambda
  round(SEM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SEM_summary$aic, 3)  # AIC
)

# SAC (SARAR)
SAC_model <- sacsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list, type = "sac")
SAC_summary <- extract_model_summary(SAC_model)

results$SAC <- c(
  format_coef(SAC_summary$coefficients[1], SAC_summary$std_errors[1], SAC_summary$p_values[1]),  # Constant
  format_coef(SAC_summary$coefficients[3], SAC_summary$std_errors[3], SAC_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(SAC_summary$coefficients[2], SAC_summary$std_errors[2], SAC_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(SAC_summary$coefficients[4], SAC_summary$std_errors[4], SAC_summary$p_values[4]),  # H_Coefficient (rel_mean)
  NA, NA, NA,  # Wx_A, Wx_E, Wx_H (Not applicable for SAC)
  round(SAC_summary$rho, 3),  # Rho
  round(SAC_summary$lambda, 3),  # Lambda
  round(SAC_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SAC_summary$aic, 3)  # AIC
)

# SLX (Spatially Lagged X Model)
SLX_model <- lmSLX(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list)
SLX_summary <- summary(SLX_model)

results$SLX <- c(
  format_coef(coef(SLX_summary)[1], SLX_summary$coefficients[1, 2], SLX_summary$coefficients[1, 4]),  # Constant
  format_coef(coef(SLX_summary)[3], SLX_summary$coefficients[3, 2], SLX_summary$coefficients[3, 4]),  # A_Coefficient (area)
  format_coef(coef(SLX_summary)[2], SLX_summary$coefficients[2, 2], SLX_summary$coefficients[2, 4]),  # E_Coefficient (elev_mean)
  format_coef(coef(SLX_summary)[4], SLX_summary$coefficients[4, 2], SLX_summary$coefficients[4, 4]),  # H_Coefficient (rel_mean)
  format_coef(coef(SLX_summary)[6], SLX_summary$coefficients[6, 2], SLX_summary$coefficients[6, 4]),  # Wx_A (lagged area)
  format_coef(coef(SLX_summary)[5], SLX_summary$coefficients[5, 2], SLX_summary$coefficients[5, 4]),  # Wx_E (lagged elev_mean)
  format_coef(coef(SLX_summary)[7], SLX_summary$coefficients[7, 2], SLX_summary$coefficients[7, 4]),  # Wx_H (lagged rel_mean)
  NA,  # Rho (Not applicable for SLX)
  NA,  # Lambda (Not applicable for SLX)
  NA,  # Nagelkerke (Not provided for SLX)
  round(AIC(SLX_model), 3)  # AIC
)

# Spatial Durbin Error Model (SDEM)
SDEM_model <- errorsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list, etype = "emixed")
SDEM_summary <- extract_model_summary(SDEM_model)

results$SDEM <- c(
  format_coef(SDEM_summary$coefficients[1], SDEM_summary$std_errors[1], SDEM_summary$p_values[1]),  # Constant
  format_coef(SDEM_summary$coefficients[3], SDEM_summary$std_errors[3], SDEM_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(SDEM_summary$coefficients[2], SDEM_summary$std_errors[2], SDEM_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(SDEM_summary$coefficients[4], SDEM_summary$std_errors[4], SDEM_summary$p_values[4]),  # H_Coefficient (rel_mean)
  format_coef(SDEM_summary$coefficients[6], SDEM_summary$std_errors[6], SDEM_summary$p_values[6]),  # Wx_A (lagged area)
  format_coef(SDEM_summary$coefficients[5], SDEM_summary$std_errors[5], SDEM_summary$p_values[5]),  # Wx_E (lagged elev_mean)
  format_coef(SDEM_summary$coefficients[7], SDEM_summary$std_errors[7], SDEM_summary$p_values[7]),  # Wx_H (lagged rel_mean)
  NA,  # Rho (Not applicable for SDEM)
  round(SDEM_summary$lambda, 3),  # Lambda
  round(SDEM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SDEM_summary$aic, 3)  # AIC
)

# Spatial Durbin Model (SDM)
SDM_model <- lagsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list, type = "mixed")
SDM_summary <- extract_model_summary(SDM_model)

results$SDM <- c(
  format_coef(SDM_summary$coefficients[1], SDM_summary$std_errors[1], SDM_summary$p_values[1]),  # Constant
  format_coef(SDM_summary$coefficients[3], SDM_summary$std_errors[3], SDM_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(SDM_summary$coefficients[2], SDM_summary$std_errors[2], SDM_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(SDM_summary$coefficients[4], SDM_summary$std_errors[4], SDM_summary$p_values[4]),  # H_Coefficient (rel_mean)
  format_coef(SDM_summary$coefficients[6], SDM_summary$std_errors[6], SDM_summary$p_values[6]),  # Wx_A (lagged area)
  format_coef(SDM_summary$coefficients[5], SDM_summary$std_errors[5], SDM_summary$p_values[5]),  # Wx_E (lagged elev_mean)
  format_coef(SDM_summary$coefficients[7], SDM_summary$std_errors[7], SDM_summary$p_values[7]),  # Wx_H (lagged rel_mean)
  round(SDM_summary$rho, 3),  # Rho
  NA,  # Lambda (Not applicable for SDM)
  round(SDM_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(SDM_summary$aic, 3)  # AIC
)

# Mansky GNS model
GNS_model <- sacsarlm(y_log ~ elev_mean + area + rel_mean, data = aoi, listw = nb_k5_list, type = "sacmixed")
GNS_summary <- extract_model_summary(GNS_model)

results$GNS <- c(
  format_coef(GNS_summary$coefficients[1], GNS_summary$std_errors[1], GNS_summary$p_values[1]),  # Constant
  format_coef(GNS_summary$coefficients[3], GNS_summary$std_errors[3], GNS_summary$p_values[3]),  # A_Coefficient (area)
  format_coef(GNS_summary$coefficients[2], GNS_summary$std_errors[2], GNS_summary$p_values[2]),  # E_Coefficient (elev_mean)
  format_coef(GNS_summary$coefficients[4], GNS_summary$std_errors[4], GNS_summary$p_values[4]),  # H_Coefficient (rel_mean)
  format_coef(GNS_summary$coefficients[6], GNS_summary$std_errors[6], GNS_summary$p_values[6]),  # Wx_A (lagged area)
  format_coef(GNS_summary$coefficients[5], GNS_summary$std_errors[5], GNS_summary$p_values[5]),  # Wx_E (lagged elev_mean)
  format_coef(GNS_summary$coefficients[7], GNS_summary$std_errors[7], GNS_summary$p_values[7]),  # Wx_H (lagged rel_mean)
  round(GNS_summary$rho, 3),  # Rho
  round(GNS_summary$lambda, 3),  # Lambda
  round(GNS_summary$nagelkerke, 3),  # Nagelkerke R^2
  round(GNS_summary$aic, 3)  # AIC
)

# Create a LaTeX table
latex_table <- xtable(results, caption = "Model Results for Various Spatial Models", label = "tab:spatial_models")

# Print LaTeX table without escaping LaTeX commands and with custom hline
print(latex_table, type = "latex", include.rownames = FALSE, hline.after = c(0, 7, 9), sanitize.text.function = identity)

####################################

# regimes-Spatial Error Model-cluster
SEM_KNN = errorsarlm(y_log ~ 0 + (elev_mean + area + rel_mean):(knn5), data = aoi, listw = nb_k5_list)
SEM_KNN_summary <- summary(SEM_KNN, Nagelkerke = TRUE)

# regimes-Spatial Lag Model-cluster
SAR_KNN = lagsarlm(y_log ~ 0 + (area + elev_mean + rel_mean):(knn5), data = aoi2, listw = nb_k5_list)
SAR_KNN_summary <- summary(SAR_KNN, Nagelkerke = TRUE)

# Extract coefficient table for both models
coeff_table_SEM <- SEM_KNN_summary$Coef
coeff_table_SAR <- SAR_KNN_summary$Coef

# Function to format coefficients with standard errors and bold significant ones
format_coef <- function(coef, std_err, p_value) {
  if (!is.na(p_value) && !is.null(p_value) && is.numeric(p_value) && p_value < 0.05) {
    return(paste0("\\textbf{", sprintf("%.3f", coef), "} (", sprintf("%.2f", std_err), ")"))
  } else {
    return(paste0(sprintf("%.3f", coef), " (", sprintf("%.2f", std_err), ")"))
  }
}

# Replace variable names with human-readable versions and clusters with full names
row_names_SEM <- rownames(coeff_table_SEM)
row_names_SEM <- gsub("elev_mean", "Mean Elevation", row_names_SEM)
row_names_SEM <- gsub("area", "Catchment Area", row_names_SEM)
row_names_SEM <- gsub("rel_mean", "Mean Relief", row_names_SEM)
row_names_SEM <- gsub("knn5A", "Cluster A", row_names_SEM)
row_names_SEM <- gsub("knn5B", "Cluster B", row_names_SEM)
row_names_SEM <- gsub("knn5C", "Cluster C", row_names_SEM)
row_names_SEM <- gsub("knn5D", "Cluster D", row_names_SEM)

# Ensure row names for SAR_KNN match the same structure for proper merging
row_names_SAR <- rownames(coeff_table_SAR)
row_names_SAR <- gsub("elev_mean", "Mean Elevation", row_names_SAR)
row_names_SAR <- gsub("area", "Catchment Area", row_names_SAR)
row_names_SAR <- gsub("rel_mean", "Mean Relief", row_names_SAR)
row_names_SAR <- gsub("knn5A", "Cluster A", row_names_SAR)
row_names_SAR <- gsub("knn5B", "Cluster B", row_names_SAR)
row_names_SAR <- gsub("knn5C", "Cluster C", row_names_SAR)
row_names_SAR <- gsub("knn5D", "Cluster D", row_names_SAR)

# Create data frames for both models with formatted coefficients
results <- data.frame(
  Covariates = row_names_SEM,
  `SEM $+$ Regimes: Cluster` = mapply(format_coef, coeff_table_SEM[, "Estimate"], coeff_table_SEM[, "Std. Error"], coeff_table_SEM[, "Pr(>|z|)"]),
  `SAR $+$ Regimes: Cluster` = mapply(format_coef, coeff_table_SAR[, "Estimate"], coeff_table_SAR[, "Std. Error"], coeff_table_SAR[, "Pr(>|z|)"])
)

# Check if lambda and lambda-related values exist in the SEM model summary
lambda_position <- NULL
if (!is.null(SEM_KNN_summary$lambda) && !is.null(SEM_KNN_summary$lambda.se)) {
  
  # Format lambda for SEM model
  lambda_SEM <- SEM_KNN_summary$lambda
  lambda_std_err_SEM <- SEM_KNN_summary$lambda.se
  lambda_p_value_SEM <- SEM_KNN_summary$lambda.pval
  
  lambda_row_SEM <- format_coef(lambda_SEM, lambda_std_err_SEM, lambda_p_value_SEM)
}

# Check if rho exists in the SAR model summary
rho_position <- NULL
if (!is.null(SAR_KNN_summary$rho) && !is.null(SAR_KNN_summary$rho.se)) {
  
  # Format rho for SAR model
  rho_SAR <- SAR_KNN_summary$rho
  rho_std_err_SAR <- SAR_KNN_summary$rho.se
  rho_p_value_SAR <- SAR_KNN_summary$rho.pval
  
  rho_row_SAR <- format_coef(rho_SAR, rho_std_err_SAR, rho_p_value_SAR)
}

# Add the lambda and rho row to the results
lambda_position <- nrow(results) + 1
results <- rbind(results, data.frame(
  Covariates = c("$\\lambda$", "$\\rho$"),
  `SEM $+$ Regimes: Cluster` = c(lambda_row_SEM, ""),  # Only SEM has lambda
  `SAR $+$ Regimes: Cluster` = c("", rho_row_SAR)  # Only SAR has rho
))

# Check if Nagelkerke R² exists, otherwise assign NaN for both models
nagelkerke_r2_SEM <- ifelse(!is.null(SEM_KNN_summary$NK), sprintf("%.4f", SEM_KNN_summary$NK), "NaN")
nagelkerke_r2_SAR <- ifelse(!is.null(SAR_KNN_summary$NK), sprintf("%.4f", SAR_KNN_summary$NK), "NaN")

# Fetch AIC values from the model objects
aic_SEM <- AIC(SEM_KNN)  # Explicitly fetch AIC from the SEM model
aic_SAR <- AIC(SAR_KNN)  # Explicitly fetch AIC from the SAR model

# Add Nagelkerke R² and AIC for both models
model_stats <- data.frame(
  Covariates = c("Nagelkerke R²", "AIC"),
  `SEM $+$ Regimes: Cluster` = c(
    nagelkerke_r2_SEM,  # SEM Nagelkerke R²
    sprintf("%.3f", aic_SEM)  # SEM AIC
  ),
  `SAR $+$ Regimes: Cluster` = c(
    nagelkerke_r2_SAR,  # SAR Nagelkerke R²
    sprintf("%.3f", aic_SAR)  # SAR AIC
  )
)

# Add model statistics to the results after Lambda and Rho
results <- rbind(results, model_stats)

# Prepare custom header using add.to.row and prevent default header from xtable
add_rows <- list(
  pos = list(0, lambda_position - 1, nrow(results) - 2, nrow(results) - 1),
  command = c(
    "\\toprule\n\\textbf{Covariates} & \\textbf{SEM $+$ Regimes: Cluster} & \\textbf{SAR $+$ Regimes: Cluster} \\\\\n\\midrule\n",
    "\\midrule\n",  # \midrule before Lambda and Rho rows
    "\\midrule\n",
    "\\bottomrule\n\\multicolumn{3}{l}{\\textbf{$p<0.05$}, (standard errors)} \\\\\n"
  )
)

# Print LaTeX table using xtable
latex_table <- xtable(results, caption = "SEM and SAR KNN Model Coefficients and Statistics", label = "tab:sem_sar_knn")

# Output the LaTeX table without default header, with the custom rules
print(latex_table, type = "latex", include.rownames = FALSE, include.colnames = FALSE,  # Prevent default header from printing
      sanitize.text.function = identity, 
      caption.placement = "top",  # Move caption to the top
      add.to.row = add_rows,
      hline.after = NULL  # Disable automatic \hline
)

