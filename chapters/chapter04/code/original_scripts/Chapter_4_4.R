# Load libraries
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gplots)) install.packages("gplots")
library(ggplot2)
library(gplots)

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

# 3. Create Heatmap visualization of Cluster Matrix
heatmap.2(as.matrix(cg_dist),
          trace = "none",
          col = colorRampPalette(c("white", "yellow", "red"))(100),
          main = "cgMLST Loci Distance Matrix",
          margins = c(8, 8),
          density.info = "none",
          cellnote = cg_dist,
          notecol = "black")