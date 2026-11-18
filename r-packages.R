# R packages used by the chapter scripts (install.packages from CRAN). Run through scripts/setup.ps1 or directly:  Rscript r-packages.R ggplot2 ape
args <- commandArgs(trailingOnly = TRUE)
pk <- if (length(args)) args else c("ggplot2", "gplots", "ape", "outbreaks", "coda", "outbreaker2", "TransPhylo")
lib <- Sys.getenv("R_LIBS_USER"); if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
need <- pk[!pk %in% rownames(installed.packages())]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org", dependencies = TRUE)
# TransPhylo is no longer on CRAN (only in the archive), so install.packages() cannot find it: fall back to the archived source (needs a C++ compiler; ape and Rcpp are installed first).
if ("TransPhylo" %in% pk && !"TransPhylo" %in% rownames(installed.packages())) {
  deps <- c("ape", "Rcpp")[!c("ape", "Rcpp") %in% rownames(installed.packages())]
  if (length(deps)) install.packages(deps, repos = "https://cloud.r-project.org")
  install.packages("https://cran.r-project.org/src/contrib/Archive/TransPhylo/TransPhylo_1.4.5.tar.gz", repos = NULL, type = "source", lib = if (nzchar(lib)) lib else .libPaths()[1])
}
bad <- pk[!pk %in% rownames(installed.packages())]
for (p in pk) cat(sprintf("%-14s %s\n", p, if (p %in% bad) "MISSING" else as.character(packageVersion(p))))
if (length(bad)) quit(status = 1)
