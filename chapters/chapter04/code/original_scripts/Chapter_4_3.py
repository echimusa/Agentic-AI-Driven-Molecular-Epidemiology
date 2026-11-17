import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# 1. Simulate core genome SNP profiles for 10 patient isolates
np.random.seed(88)
n_isolates = 10
n_snp_loci = 50  # Variable positions in core alignment

# Create labels (e.g., patient samples from a suspected outbreak)
isolate_labels = [f"Pat_{i:02d}" for i in range(1, n_isolates + 1)]

# Establish an ancestral core genome sequence (0 = ancestral allele)
ancestral = np.zeros(n_snp_loci)

# Define sub-cluster groups to simulate transmission lineages
# Group 1 (Outbreak Cluster A - active transmission): 4 closely related isolates
# Group 2 (Outbreak Cluster B - separate introduction): 3 closely related isolates
# Diverse background cases (unrelated endemic isolates): 3 isolated strains
genotypes = np.zeros((n_isolates, n_snp_loci))

# Lineage A (Simulated outbreak, tightly clustered)
genotypes[0] = ancestral + np.random.binomial(1, 0.02, n_snp_loci) # Patient 1
for i in range(1, 4):
    # Introduce 0-2 micro-evolutionary SNPs from Case 1 (direct transmission)
    genotypes[i] = genotypes[0] ^ np.random.binomial(1, 0.04, n_snp_loci)

# Lineage B (Secondary cluster, distinct lineage)
genotypes[4] = ancestral + np.random.binomial(1, 0.35, n_snp_loci) # Patient 5
for i in range(5, 7):
    # Clustered together, far from Ancestral and Lineage A
    genotypes[i] = genotypes[4] ^ np.random.binomial(1, 0.04, n_snp_loci)

# Unrelated background endemic cases
genotypes[7] = np.random.binomial(1, 0.40, n_snp_loci)
genotypes[8] = np.random.binomial(1, 0.45, n_snp_loci)
genotypes[9] = np.random.binomial(1, 0.50, n_snp_loci)

# 2. Calculate Pairwise SNP Distance Matrix
snp_matrix = np.zeros((n_isolates, n_isolates), dtype=int)
for i in range(n_isolates):
    for j in range(n_isolates):
        # Manhattan distance on binary arrays calculates total base differences
        snp_matrix[i, j] = int(np.sum(np.abs(genotypes[i] - genotypes[j])))

snp_df = pd.DataFrame(snp_matrix, index=isolate_labels, columns=isolate_labels)

print("--- Pairwise SNP Distance Matrix (Delta SNPs) ---")
print(snp_df)

# Evaluate Genomic Outbreak threshold (typically <= 5-12 SNPs for bacterial pathogens)
print("\nIsolated Pairs Passing Outbreak Closeness Threshold (<= 12 SNPs):")
pairs_found = []
for i in range(n_isolates):
    for j in range(i+1, n_isolates):
        if snp_matrix[i, j] <= 12:
            print(f"Possible Transmission Match: {isolate_labels[i]} <-> {isolate_labels[j]} ({snp_matrix[i,j]} SNPs)")
            pairs_found.append((isolate_labels[i], isolate_labels[j], snp_matrix[i,j]))

# 3. Perform Hierarchical Clustering and Visualization
condensed_dist = squareform(snp_matrix)
linkage_matrix = linkage(condensed_dist, method='average')

plt.figure(figsize=(10, 5))
# Generate Clustered Heatmap
sns.heatmap(snp_df, annot=True, cmap="YlOrRd", cbar_kws={'label': 'SNP Differences'})
plt.title("SNP Pairwise Distance and Clustering Matrix")
plt.tight_layout()
plt.show()

# Visualize Clustering Trees (Dendrograms)
plt.figure(figsize=(7, 4))
dendrogram(linkage_matrix, labels=isolate_labels, color_threshold=12)
plt.axhline(y=12, color='red', linestyle='--', label='Outbreak Cutoff (12 SNPs)')
plt.title("Isolate Clustering via Hierarchical Tree")
plt.ylabel("Genetic Distance (SNP Count)")
plt.legend()
plt.tight_layout()
plt.show()