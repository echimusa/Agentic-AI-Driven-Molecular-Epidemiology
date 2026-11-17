# Chapter 1, Hands-on 2 (R check): pairwise SNP distances and their correlation with sampling time, in base R.
# Inputs: ebola_aln.fasta (aligned genomes) and ebola_dates.csv (accession, sampling date).
lines <- readLines("ebola_aln.fasta"); hdr <- which(startsWith(lines, ">")); ids <- sub("^>(\\S+).*", "\\1", lines[hdr]); ends <- c(hdr[-1] - 1, length(lines))
seqs <- mapply(function(a, b) toupper(paste(lines[(a + 1):b], collapse = "")), hdr, ends); M <- do.call(rbind, strsplit(seqs, "")); valid <- M %in% c("A", "C", "G", "T"); dim(valid) <- dim(M)
n <- nrow(M); D <- matrix(0L, n, n)
for (i in seq_len(n)) D[i, ] <- rowSums((M != matrix(M[i, ], n, ncol(M), byrow = TRUE)) & valid & matrix(valid[i, ], n, ncol(M), byrow = TRUE))
dates <- read.csv("ebola_dates.csv"); day <- as.numeric(as.Date(dates$date[match(ids, dates$acc)])); up <- upper.tri(D)
snps <- D[up]; gap <- abs(outer(day, day, "-"))[up]
cat(sprintf("%d genomes; %d pairs; median distance %d SNPs; pairs within 2 SNPs %.1f%%; Spearman rho %.3f\n", n, length(snps), as.integer(median(snps)), 100 * mean(snps <= 2), cor(gap, snps, method = "spearman")))
