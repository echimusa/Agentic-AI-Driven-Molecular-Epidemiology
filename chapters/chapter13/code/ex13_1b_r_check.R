# Chapter 13, Hands-on 1 (R check): the growth rate and one Cori reproduction number recomputed in base R.
# Input: ebola_sl_linelist.csv (exported from the CRAN package outbreaks). Needs base R only.
d <- read.csv("ebola_sl_linelist.csv"); d$date <- as.Date(d$date_of_onset)
day <- as.data.frame(table(factor(d$date, levels = as.character(seq(min(d$date), max(d$date), by = "day")))), stringsAsFactors = FALSE); names(day) <- c("date", "n"); day$date <- as.Date(day$date)
w_gamma <- function(mean, sd, n = 60) { k <- (mean / sd)^2; th <- sd^2 / mean; e <- 0.5 + 0:n; p <- diff(pgamma(c(0, e[-1]), shape = k, scale = th)); p / sum(p) }
w <- w_gamma(15.3, 9.3)
# growth rate, 1 June to 15 September 2014: quasi-Poisson log-linear regression
s <- subset(day, date >= as.Date("2014-06-01") & date <= as.Date("2014-09-15")); s$t <- seq_len(nrow(s)); fit <- glm(n ~ t, family = quasipoisson, data = s)
r <- coef(fit)["t"]; se <- summary(fit)$coefficients["t", "Std. Error"]; k <- (15.3 / 9.3)^2; th <- 9.3^2 / 15.3
cat(sprintf("growth rate r = %.4f per day (SE %.4f); doubling %.1f d; R = %.2f\n", r, se, log(2) / r, (1 + r * th)^k))
# Cori (2013) posterior mean of Rt over the 14 days ending 1 October 2014, prior Gamma(1, scale 5)
inc <- day$n; tt <- which(day$date == as.Date("2014-10-01")); lam <- sapply(seq_along(inc), function(t) if (t == 1) 0 else { m <- min(t - 1, length(w)); sum(w[1:m] * inc[(t - 1):(t - m)]) })
win <- (tt - 13):tt; shape <- 1 + sum(inc[win]); rate <- 1 / 5 + sum(lam[win]); cat(sprintf("Rt on 1 Oct 2014 = %.2f (95%% interval %.2f to %.2f)\n", shape / rate, qgamma(0.025, shape, rate), qgamma(0.975, shape, rate)))
