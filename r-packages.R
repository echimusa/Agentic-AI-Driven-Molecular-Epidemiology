# R packages used by the chapter scripts (install.packages from CRAN). Run through scripts/setup.ps1 or directly:  Rscript r-packages.R ggplot2 ape
args <- commandArgs(trailingOnly = TRUE)
pk <- if (length(args)) args else c("ggplot2", "gplots", "ape", "outbreaks", "coda", "outbreaker2", "TransPhylo")
lib <- Sys.getenv("R_LIBS_USER"); if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
need <- pk[!pk %in% rownames(installed.packages())]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org", dependencies = TRUE)
bad <- pk[!pk %in% rownames(installed.packages())]
for (p in pk) cat(sprintf("%-14s %s\n", p, if (p %in% bad) "MISSING" else as.character(packageVersion(p))))
if (length(bad)) quit(status = 1)
