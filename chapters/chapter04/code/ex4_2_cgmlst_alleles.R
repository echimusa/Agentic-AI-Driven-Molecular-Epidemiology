# Load libraries
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

# 1. Generate Simulated cgMLST Allelic Profiles for 12 Isolates (20 Core Loci)
set.seed(42)
n_isolates <- 12
n_loci <- 20 # 20 Core gene loci tracking allelic shifts

isolate_names <- paste0("Isolate_", 1:n_isolates)
locus_names <- paste0("Locus_", 1:n_loci)

# Define baseline allele values (using standard arbitrary integer identifiers)
alleles <- matrix(sample(1:5, n_isolates * n_loci, replace = TRUE), 
                  nrow = n_isolates, ncol = n_loci)

# Introduce Cluster 1 (Active Hospital Outbreak: 0-1 allelic locus shifts)
# Isolate_1, 2, 3, 4 represent a localized outbreak
for (i in 2:4) {
  alleles[i, ] <- alleles[1, ] # start identical to Isolate_1
  # Randomly alter 1 locus to show natural shift
  mutated_locus <- sample(1:n_loci, 1)
  alleles[i, mutated_locus] <- sample(6:10, 1) # Assign brand new allele
}

# Introduce Cluster 2 (Community transmission: closely grouped independently)
# Isolate_5, 6, 7
for (i in 6:7) {
  alleles[i, ] <- alleles[5, ]
  mutated_loci <- sample(1:n_loci, 2)
  alleles[i, mutated_loci] <- sample(11:15, 2) # assign distinct alleles
}

cg_profiles <- as.data.frame(alleles)
rownames(cg_profiles) <- isolate_names
colnames(cg_profiles) <- locus_names

print("--- Simulated cgMLST Profiles (Individual Locus Allele Numbers) ---")
print(head(cg_profiles, 5))

# 2. Compute Pairwise cgMLST Allelic Distances (Mismatch Fraction * n_loci)
# Count exactly how many core genes carry different alleles
calc_cg_distance <- function(profile_df) {
  n <- nrow(profile_df)
  dist_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist_matrix[i, j] <- sum(profile_df[i, ] != profile_df[j, ])
    }
  }
  rownames(dist_matrix) <- rownames(profile_df)
  colnames(dist_matrix) <- rownames(profile_df)
  return(dist_matrix)
}

cg_dist <- calc_cg_distance(cg_profiles)

print("\n--- cgMLST Pairwise Distance Matrix (Locus Mismatches) ---")
print(cg_dist[1:6, 1:6])

# Set a typical strict transmission cutoff: matching alleles at >= 85% of loci
# with 20 loci, this allows <= 3 allele differences
cutoff <- 3
outbreaks <- which(cg_dist <= cutoff & lower.tri(cg_dist), arr.ind = TRUE)

print(paste("\nActive Clusters under Cluster Cutoff (<=", cutoff, "alleles):"))
for(idx in 1:nrow(outbreaks)) {
  i <- outbreaks[idx, 1]
  j <- outbreaks[idx, 2]
  print(paste("Epidemiological Link Identified:", 
              isolate_names[i], "<->", isolate_names[j], 
              "differing by only", cg_dist[i, j], "alleles"))
}

# 3. Additions for the chapter: pairs against clusters. A pairwise cut-off links two isolates only when they are close to each other; a clustering rule also joins isolates through a third one.
single <- hclust(as.dist(cg_dist), method = "single"); complete <- hclust(as.dist(cg_dist), method = "complete"); average <- hclust(as.dist(cg_dist), method = "average")
cl_single <- cutree(single, h = cutoff); cl_complete <- cutree(complete, h = cutoff); cl_average <- cutree(average, h = cutoff)
cat("\nClusters at a cut-off of", cutoff, "alleles\n"); print(data.frame(isolate = isolate_names, single = cl_single, average = cl_average, complete = cl_complete), row.names = FALSE)
cat(sprintf("\nNumber of clusters with 2 or more isolates: single %d, average %d, complete %d\n", sum(table(cl_single) > 1), sum(table(cl_average) > 1), sum(table(cl_complete) > 1)))
same <- outer(cl_single, cl_single, "==") & lower.tri(cg_dist); chained <- which(same & cg_dist > cutoff, arr.ind = TRUE)
cat("Pairs in the same single-linkage cluster although more than", cutoff, "alleles apart:", nrow(chained), "\n")
if (nrow(chained) > 0) for (k in seq_len(nrow(chained))) cat(sprintf("  %s and %s differ by %d alleles\n", isolate_names[chained[k, 1]], isolate_names[chained[k, 2]], cg_dist[chained[k, 1], chained[k, 2]]))
cat(sprintf("Mean allele distance within the simulated outbreaks: %.1f; between unrelated isolates: %.1f (of %d loci)\n", mean(cg_dist[1:4, 1:4][lower.tri(cg_dist[1:4, 1:4])]), mean(cg_dist[8:12, 8:12][lower.tri(cg_dist[8:12, 8:12])]), n_loci))
# 4. Heatmap with the base R function (the uploaded script used gplots::heatmap.2; this avoids the extra package)
dir.create("figs4", showWarnings = FALSE); png("figs4/fig2_cgmlst_heatmap_R.png", width = 1500, height = 1300, res = 200)
heatmap(as.matrix(cg_dist), Rowv = as.dendrogram(average), Colv = as.dendrogram(average), col = colorRampPalette(c("white", "yellow", "red"))(100), scale = "none", margins = c(8, 8), main = "cgMLST allele distance matrix (12 simulated isolates)")
dev.off()
