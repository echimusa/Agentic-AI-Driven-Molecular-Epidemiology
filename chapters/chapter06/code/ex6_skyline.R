# Chapter 6, Exercise 2 - A coalescent skyline plot, checked against a known truth.
# We simulate a genealogy under exponential population growth, estimate the population size
# from the coalescent intervals, and compare with the size we simulated. Base R + ggplot2 only.
library(ggplot2)
set.seed(42)

N0 <- 1000      # effective population size at the present (in units of "N x generation time", years)
r  <- 0.5       # growth rate per year: N(t) = N0 * exp(-r t), with t in years BEFORE the present
n  <- 100       # number of sampled lineages

# Simulate coalescent times backwards in time under exponential growth (inversion method).
# With k lineages the coalescent rate is k(k-1) / (2 N(t)).
sim_coalescent <- function(n, N0, r) {
  t <- 0; times <- numeric(n - 1)
  for (k in n:2) {
    rate0 <- k * (k - 1) / (2 * N0)
    t <- log(exp(r * t) - r * log(runif(1)) / rate0) / r
    times[n - k + 1] <- t
  }
  times                                              # times of the n-1 coalescent events
}
ct <- sim_coalescent(n, N0, r)

# Classic skyline estimator (Pybus et al. 2000): in the interval with k lineages of length g_k,
# N_hat = g_k * k(k-1) / 2. The generalised skyline pools m consecutive intervals to reduce noise.
k  <- n:2
g  <- diff(c(0, ct))
sky <- data.frame(t = c(0, head(ct, -1)) + g / 2, N = g * k * (k - 1) / 2, kind = "classic skyline (every interval)")
m  <- 10
grp <- ceiling(seq_along(k) / m)
gen <- do.call(rbind, lapply(split(seq_along(k), grp), function(i)
  data.frame(t = mean(c(0, ct)[c(min(i), max(i) + 1)]), N = sum(g[i] * k[i] * (k[i] - 1) / 2) / length(i),
             kind = "generalised skyline (10 intervals per group)")))
est <- rbind(sky, gen)
tt <- seq(0, max(ct), length.out = 200)
truth <- data.frame(t = tt, N = N0 * exp(-r * tt))

# Recover the growth rate from the generalised skyline by regressing log N on t
fit <- lm(log(N) ~ t, data = gen)
cat(sprintf("simulated: N0 = %d, r = %.2f per year, %d lineages, root at %.2f years\n", N0, r, n, max(ct)))
cat(sprintf("recovered growth rate from the generalised skyline: %.2f per year (95%% CI %.2f to %.2f)\n",
            -coef(fit)[["t"]], -confint(fit)["t", 2], -confint(fit)["t", 1]))
cat(sprintf("recovered current size (intercept): %.0f (true %d)\n", exp(coef(fit)[["(Intercept)"]]), N0))

# Parametric alternative: maximum likelihood for an exponential-growth coalescent (what the exponential
# tree prior in BEAST 2 fits). log L = sum[ log(k(k-1)/2) - log N(t_i) ] - sum[ k(k-1)/2 * integral of 1/N ]
mle_exp <- function(ct, n) {
  k <- n:2; tprev <- c(0, head(ct, -1)); kk <- k * (k - 1) / 2
  nll <- function(par) {                              # par = (log N0, log r)
    lN0 <- par[1]; rr <- exp(par[2])
    -(sum(log(kk) + rr * ct - lN0) - sum(kk * (exp(rr * ct) - exp(rr * tprev)) / (rr * exp(lN0))))
  }
  o <- optim(c(log(500), log(0.3)), nll, method = "Nelder-Mead")
  c(N0 = exp(o$par[1]), r = exp(o$par[2]))
}
mle <- mle_exp(ct, n)
cat(sprintf("maximum-likelihood exponential fit on the same genealogy: N0 = %.0f, r = %.2f per year\n", mle[["N0"]], mle[["r"]]))

# Repeat the whole experiment 200 times to see bias and spread of the two estimators
rep <- t(replicate(200, {
  c2 <- sim_coalescent(n, N0, r); g2 <- diff(c(0, c2)); k2 <- n:2
  gr <- ceiling(seq_along(k2) / m)
  tt2 <- tapply(seq_along(k2), gr, function(i) mean(c(0, c2)[c(min(i), max(i) + 1)]))
  nn2 <- tapply(seq_along(k2), gr, function(i) sum(g2[i] * k2[i] * (k2[i] - 1) / 2) / length(i))
  c(skyline = -coef(lm(log(nn2) ~ tt2))[[2]], mle = mle_exp(c2, n)[["r"]])
}))
cat(sprintf("over 200 simulated genealogies (true r = %.2f): skyline regression %.2f (SD %.2f); exponential ML %.2f (SD %.2f)\n",
            r, mean(rep[, "skyline"]), sd(rep[, "skyline"]), mean(rep[, "mle"]), sd(rep[, "mle"])))

p <- ggplot() +
  geom_point(data = subset(est, kind == "classic skyline (every interval)"), aes(t, N, colour = kind), alpha = 0.45, size = 1.6) +
  geom_step(data = subset(est, kind != "classic skyline (every interval)"), aes(t, N, colour = kind), linewidth = 1) +
  geom_line(data = truth, aes(t, N, colour = "true population size"), linewidth = 1.1, linetype = 2) +
  scale_y_log10() +
  scale_colour_manual(values = c("classic skyline (every interval)" = "#718096",
                                 "generalised skyline (10 intervals per group)" = "#2b6cb0",
                                 "true population size" = "#c05621"), name = NULL) +
  labs(x = "Years before the present", y = "Effective population size (log scale)",
       title = "Skyline estimates track exponential growth from 100 simulated lineages") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
dir.create("figs6", showWarnings = FALSE)
ggsave("figs6/fig_ch6_skyline.png", p, width = 8, height = 4.6, dpi = 200)
cat("figure written to figs6/fig_ch6_skyline.png\n")
