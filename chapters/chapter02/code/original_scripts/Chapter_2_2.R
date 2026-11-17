# Load libraries (install first if needed)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

# 1. Simulate Biomarker Data
set.seed(123)
n_samples <- 100
n_markers <- 20 # 20 molecular assays (e.g., proteins or transcripts)

# Simulate natural biological variation
base_matrix <- matrix(rnorm(n_samples * n_markers, mean = 5, sd = 1), 
                      nrow = n_samples, ncol = n_markers)

# Assign samples to laboratory Batch 1 or Batch 2
batches <- factor(c(rep("Batch_1", 50), rep("Batch_2", 50)))

# Artificially inject a shift into Batch 2 (Batch Effect)
batch_effect <- rep(0, n_samples)
batch_effect[batches == "Batch_2"] <- 1.5 # systematic artifact shift

# Add batch bias to all molecular assays
for(j in 1:n_markers) {
  base_matrix[, j] <- base_matrix[, j] + batch_effect
}

biomarker_df <- as.data.frame(base_matrix)
colnames(biomarker_df) <- paste0("Marker_", 1:n_markers)

# 2. Run PCA to detect Batch Effect
pca_results <- prcomp(biomarker_df, scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  Batch = batches
)

# Render Scatter Plot to visualize the bias
ggplot(pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "PCA Plot for Pre-Analytical Batch Effect Detection",
    subtitle = "Distinct clustering shows artificial batch variation",
    x = paste0("PC1 (Variance: ", round(100 * summary(pca_results)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (Variance: ", round(100 * summary(pca_results)$importance[2, 2], 1), "%)")
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e"))