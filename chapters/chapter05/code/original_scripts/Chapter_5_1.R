library(ape)
library(phangorn)
library(ggplot2)

# 1. Download and read 15 Influenza H3N2 sequences (FASTA) from a suspected outbreak
# This dataset contains aligned HA gene sequences.
fasta_url <- "https://raw.githubusercontent.com/colinmeir/molecular_epi/main/data/influenza_h3n2_outbreak.fasta"
temp_file <- tempfile(fileext = ".fasta")
download.file(fasta_url, temp_file, quiet = TRUE)

# Load the alignment
dna_align <- read.dna(temp_file, format = "fasta")
print(paste("Dataset loaded:", nrow(dna_align), "sequences of length", ncol(dna_align), "bp"))

# 2. Compute Genetic Distance using Jukes-Cantor model (handling observed transitions/transversions)
jc_distance <- dist.dna(dna_align, model = "JC69")
print("Preview of pairwise genetic distance matrix:")
print(round(as.matrix(jc_distance)[1:5, 1:5], 5))

# 3. Construct Neighbor-Joining Tree
nj_tree <- nj(jc_distance)

# 4. Assess Node Support via Bootstrap Resampling (100 replicates)
# Bootstrap evaluates the reliability of tree branching patterns
set.seed(42)
boot_trees <- boot.phylo(nj_tree, dna_align, function(x) nj(dist.dna(x, model = "JC69")), B = 100)

# Map bootstrap scores back to tree nodes
nj_tree$node.label <- round(boot_trees, 0)

# Plot the labeled phylogenetic tree
plot(nj_tree, 
     main = "Reconstructed Influenza A H3N2 Outbreak Phylogeny", 
     sub = "NJ Tree with Bootstrap Values (100 Reps)", 
     cex = 0.8, 
     edge.width = 1.5,
     no.margin = FALSE)

# Add bootstrap values to branching nodes (support scores > 70 are typically considered reliable)
nodelabels(nj_tree$node.label, adj = c(1.2, -0.3), frame = "none", bg = "transparent", col = "blue", cex = 0.7)
add.scale.bar(length = 0.005) # branch scale represents estimated mutations per site