# Chapter 7, R companion 2 - TransPhylo on simulated outbreaks with a known truth.
# For each sampling probability pi, simulate outbreaks with the Session 7 parameters, infer the transmission tree from the sampled phylogeny alone,
# and write what Python needs to score it against the truth: <prefix>_runs.csv, _truth.csv, _pairs.csv. Environment variables TP_ITER, TP_THIN and TP_PREFIX set the chain length, the thinning and the file prefix (defaults 100000, 10, transphylo; the short-chain comparison uses 2000, 1, transphylo_short). Needs TransPhylo 1.4.5 (archived on CRAN; install from the CRAN archive) and coda.
suppressMessages({library(TransPhylo); library(coda)})
PIS <- c(0.25, 0.75); REPS <- 6; ITER <- as.integer(Sys.getenv("TP_ITER", "100000")); THIN <- as.integer(Sys.getenv("TP_THIN", "10")); PREFIX <- Sys.getenv("TP_PREFIX", "transphylo"); NMIN <- 15; NMAX <- 35
NEG <- 100/365; OFFR <- 5; WSH <- 10; WSC <- 0.1; START <- 2005; DATET <- 2008
sim_one <- function(pi, seed) {
  for (a in 0:500) {
    set.seed(seed * 1000 + a)
    sim <- suppressMessages(simulateOutbreak(neg = NEG, pi = pi, off.r = OFFR, off.p = 0.5, w.shape = WSH, w.scale = WSC, dateStartOutbreak = START, dateT = DATET))
    tt <- extractTTree(sim)$ttree; ns <- sum(!is.na(tt[, 2]))
    if (ns >= NMIN && ns <= NMAX) return(list(sim = sim, tt = tt, ns = ns, attempt = a))
  }
  stop("no outbreak of the wanted size")
}
runs <- list(); truth <- list(); pairs <- list(); k <- 0
for (pi in PIS) for (r in 1:REPS) {
  k <- k + 1; f <- sprintf("%s_run%02d.rds", PREFIX, k)
  if (file.exists(f)) { z <- readRDS(f); runs[[k]] <- z$runs; truth[[k]] <- z$truth; pairs[[k]] <- z$pairs; next }      # a finished run is not repeated after a restart
  o <- sim_one(pi, k); ptree <- extractPTree(o$sim); t0 <- Sys.time()
  res <- inferTTree(ptree, mcmcIterations = ITER, thinning = THIN, w.shape = WSH, w.scale = WSC, dateT = DATET, verbose = FALSE)
  n <- length(res); b <- floor(n / 2); post <- (b + 1):n
  lp <- sapply(res[post], function(x) x$pTTree)
  runs[[k]] <- data.frame(pi_true = pi, rep = r, run = k, seed = k, attempt = o$attempt, hosts = nrow(o$tt), sampled = o$ns, iterations = ITER, thinning = THIN, samples = n, burnin_samples = b,
    ess_log_pTTree = as.numeric(effectiveSize(lp)), pi_post = mean(sapply(res[post], function(x) x$pi)), off_r_post = mean(sapply(res[post], function(x) x$off.r)), neg_post = mean(sapply(res[post], function(x) x$neg)), seconds = as.numeric(Sys.time() - t0, units = "secs"))
  truth[[k]] <- data.frame(run = k, host = seq_len(nrow(o$tt)), infection_time = o$tt[, 1], sampling_time = o$tt[, 2], infector = o$tt[, 3])
  m <- computeMatWIW(res, burnin = 0.5); ii <- rep(seq_len(o$ns), times = o$ns); jj <- rep(seq_len(o$ns), each = o$ns)
  pairs[[k]] <- data.frame(run = k, i = ii, j = jj, prob = as.vector(m))
  saveRDS(list(runs = runs[[k]], truth = truth[[k]], pairs = pairs[[k]]), f)
  cat(sprintf("run %d (pi %.2f, rep %d): %d hosts, %d sampled, ESS %.0f, %.0f s\n", k, pi, r, nrow(o$tt), o$ns, runs[[k]]$ess_log_pTTree, runs[[k]]$seconds)); flush.console()
}
write.csv(do.call(rbind, runs), paste0(PREFIX, "_runs.csv"), row.names = FALSE); write.csv(do.call(rbind, truth), paste0(PREFIX, "_truth.csv"), row.names = FALSE)
p <- do.call(rbind, pairs); p$prob <- round(p$prob, 5); write.csv(p, paste0(PREFIX, "_pairs.csv"), row.names = FALSE); cat("written\n")
