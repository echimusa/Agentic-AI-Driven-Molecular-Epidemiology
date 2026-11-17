# Chapter 9, Hands-on 5 (extension) - the live demo of Session 9 in R: load the AMR calls, count them by drug class, compare with RGI.
# Input: amrfinder_outputs/*.amrfinder.tsv (ex9_5) and rgi_outputs/*.rgi.txt (RGI 6.0.8, CARD 4.0.2).  Output: table_class_counts_r.csv
files <- list.files("amrfinder_outputs", pattern = "amrfinder.tsv$", full.names = TRUE)
amr <- do.call(rbind, lapply(files, function(f) { d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE); d$assembly <- rep(substr(basename(f), 1, 15), nrow(d)); d }))
core <- subset(amr, Type == "AMR" & Scope == "core" & !grepl("^PARTIAL", Method) & !grepl("^INTERNAL_STOP", Method))   # the elements kept in ex9_5
cat("AMRFinderPlus core rows:", nrow(core), "\n")                                                             # an element found on two contigs is listed twice
core <- core[!duplicated(core[, c("assembly", "Element symbol")]), ]                                              # one row per element and genome, as in ex9_5
cat("distinct elements:", nrow(core), "in", length(unique(core$assembly)), "of", length(files), "genomes\n")
print(table(core$Class)); print(table(core$Subclass))
rf <- list.files("rgi_outputs", pattern = "rgi.txt$", full.names = TRUE); rf <- rf[!grepl("loose", rf)]            # the second run with Loose hits is used in ex9_5b only
rgi <- do.call(rbind, lapply(rf, function(f) { d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE); d$assembly <- rep(substr(basename(f), 1, 15), nrow(d)); d }))
cat("\nRGI hits:", nrow(rgi), " Perfect:", sum(rgi$Cut_Off == "Perfect"), " Strict:", sum(rgi$Cut_Off == "Strict"), "\n")
per <- merge(as.data.frame(table(assembly = core$assembly), responseName = "amrfinder"), as.data.frame(table(assembly = rgi$assembly), responseName = "rgi"), all = TRUE)
per$amrfinder[is.na(per$amrfinder)] <- 0; print(per)           # a genome with no core element is absent from the first table and counts as 0
cls <- as.data.frame(table(Class = core$Class), responseName = "elements"); write.csv(cls, "table_class_counts_r.csv", row.names = FALSE)
