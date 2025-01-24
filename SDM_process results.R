library(data.table)
library(dplyr)
library(terra)
library(spocc)
library(sf)
library(ggplot2)
library(readr)
library(rnaturalearth)
library(ggplot2)
library(data.table)
library(terra)
library(ggplot2)
library(SSDM)

setwd("D:/PhD/SDM/Arun")

# Load Lantana camara Model
load("Results/Lantana_camara.rda")
model <- get( "ESDM_Prosopis_juliflora")
output_dir <- "Results/Outputs/"
#plot(ESDM_Lantana_camara)
# List All Available Slots in the Model
all_slots <- slotNames(model)
print(all_slots)


#results----
##1. variable importance----
evaluation_metrics <- ESDM_Lantana_camara@evaluation
print(evaluation_metrics)
write.csv(as.data.frame(evaluation_metrics), "Results/evaluation_metrics.csv")
model_predictions <- ESDM_Lantana_camara@projection
plot(model_predictions)
dev.copy(png, filename = "model_predictions.png")
dev.off()

variable_importance <- as.data.frame(ESDM_Lantana_camara@variable.importance)

if(ncol(variable_importance) > 1) {
  variable_importance <- data.frame(
    Variable = rep(names(variable_importance), each = nrow(variable_importance)),
    Importance = unlist(variable_importance, use.names = FALSE)
  )
}

new_variable_names <- c(
  "Annual Mean Temperature",
  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
  "Isothermality (BIO2/BIO7) (* 100)",
  "Temperature Seasonality (standard deviation *100)",
  "Max Temperature of Warmest Month",
  "Min Temperature of Coldest Month",
  "Temperature Annual Range (BIO5-BIO6)",
  "Mean Temperature of Wettest Quarter",
  "Mean Temperature of Driest Quarter",
  "Mean Temperature of Warmest Quarter",
  "Mean Temperature of Coldest Quarter",
  "Annual Precipitation",
  "Precipitation of Wettest Month",
  "Precipitation of Driest Month",
  "Precipitation Seasonality (Coefficient of Variation)",
  "Precipitation of Wettest Quarter",
  "Precipitation of Driest Quarter",
  "Precipitation of Warmest Quarter",
  "Precipitation of Coldest Quarter",
  "Elevation"
)
if(length(new_variable_names) == nrow(variable_importance)) {
  variable_importance$Variable <- new_variable_names
} else {
  stop("The number of new variable names does not match the number of rows in the variable_importance dataframe.")
}

ggplot(variable_importance, aes(x = reorder(Variable, Importance), y = Importance, fill = Variable)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEEAD", 
    "#D4A5A5", "#9B59B6", "#3498DB", "#E74C3C", "#2ECC71",
    "#F1C40F", "#E67E22", "#16A085", "#8E44AD", "#2980B9",
    "#C0392B", "#27AE60", "#D35400", "#7F8C8D", "#2C3E50"
  )) +
  labs(title = "Variable Importance Lantana camara L.", x = "Variables", y = "Importance") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("Results/variable_importance_LC.png")

##2. Between-Algorithm Correlation----
algorithm_correlation <- model@algorithm.correlation
write.csv(as.data.frame(algorithm_correlation), paste0(output_dir, "Prosopis_juliflora_algorithm_correlation.csv"))
print("Between-algorithm correlation saved.")
# Visualization for Algorithm Correlation (Heatmap)
# Ensure algorithm_correlation is a matrix and includes proper row and column names
if (!is.matrix(algorithm_correlation)) {
  algorithm_correlation <- as.matrix(algorithm_correlation)
}

# Set row and column names to match the algorithms
algorithm_names <- c('SVM', 'MARS', 'RF', 'GLM', 'ANN', 'CTA', 'GBM')
rownames(algorithm_correlation) <- algorithm_names
colnames(algorithm_correlation) <- algorithm_names

# Reshape the matrix into long format
library(reshape2)
corr_df <- melt(algorithm_correlation)

# Rename columns for clarity
colnames(corr_df) <- c("Algorithm1", "Algorithm2", "Correlation")

# Plot the heatmap with AUC values
library(ggplot2)
ggplot(corr_df, aes(x = Algorithm1, y = Algorithm2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 4) + # Add AUC values
  labs(
    title = "Algorithm Correlation with AUC Values - Prosopis juliflora (Sw.) DC.",
    x = "Algorithm", y = "Algorithm", fill = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
ggsave(paste0(output_dir, "Prosopis_juliflora_algorithm_correlation_heatmap.png"))

##3. Evaluation Metrics for Thresholding----
threshold_metrics <- ESDM_Lantana_camara@algorithm.evaluation
write.csv(as.data.frame(threshold_metrics), paste0(output_dir, "Lantana_camara_threshold_metrics.csv"))
print("Evaluation metrics for thresholding saved.")

##4. Binary Map----
binary_map <- model@binary
plot(binary_map, main = "Binary Map - Prosopis juliflora (Sw.) DC.")
writeRaster(binary_map, paste0(output_dir, "Prosopis_juliflora_binary_map.tif"), format = "GTiff", overwrite = TRUE)
dev.copy(png, filename = paste0(output_dir, "Prosopis_juliflora_binary_map.png"))
dev.off()
