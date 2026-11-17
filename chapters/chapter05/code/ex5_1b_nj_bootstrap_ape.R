# Chapter 5, Hands-on 1 (R): the uploaded Chapter_5_1.R rewritten for the real H3N2 alignment of this chapter.
# Changes: the alignment is read from the local file (the uploaded URL returns HTTP 404); phangorn is not needed and is dropped; the bootstrap values are drawn on the tree with ape only.
library(ape)
dna <- read.dna("h3n2_ha_aln.fasta", format = "fasta"); cat("Dataset loaded:", nrow(dna), "sequences of length", ncol(dna), "positions\n")
jc <- dist.dna(dna, model = "JC69", pairwise.deletion = TRUE); k80 <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE); p <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)
cat(sprintf("largest distance: p %.4f, JC69 %.4f, K80 %.4f\n", max(p), max(jc), max(k80)))
nj_tree <- nj(jc); set.seed(42)
bs <- boot.phylo(nj_tree, dna, function(x) nj(dist.dna(x, model = "JC69", pairwise.deletion = TRUE)), B = 100)
cat(sprintf("NJ bootstrap (100 replicates) on the %d internal nodes: >= 95: %.1f%%, >= 70: %.1f%%, < 50: %.1f%%\n", sum(!is.na(bs)), 100 * mean(bs >= 95, na.rm = TRUE), 100 * mean(bs >= 70, na.rm = TRUE), 100 * mean(bs < 50, na.rm = TRUE)))
dir.create("figs5", showWarnings = FALSE); png("figs5/fig1b_nj_bootstrap_R.png", width = 1500, height = 1700, res = 200)
yr <- as.integer(sub(".*,", "", sub("\"", "", readLines("h3n2_meta.csv")[-1]))); nm <- sub(",.*", "", readLines("h3n2_meta.csv")[-1]); cols <- c("#2b6cb0", "#2f855a", "#d69e2e", "#dd6b20", "#c53030")[yr[match(nj_tree$tip.label, nm)] - 2013]
plot(nj_tree, cex = 0.55, tip.color = cols, edge.width = 1.2, main = "H3N2 HA: NJ tree (JC69) with bootstrap values", no.margin = FALSE); nodelabels(ifelse(!is.na(bs) & bs >= 70, bs, ""), frame = "none", col = "blue", cex = 0.55, adj = c(1.2, -0.3)); add.scale.bar(length = 0.005)
legend("bottomleft", legend = 2014:2018, col = c("#2b6cb0", "#2f855a", "#d69e2e", "#dd6b20", "#c53030"), pch = 16, cex = 0.7, bty = "n"); invisible(dev.off())
