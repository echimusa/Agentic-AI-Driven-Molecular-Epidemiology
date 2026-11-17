# Chapter 11, Hands-on 2 (R companion): the IVW, MR-Egger and weighted-median estimates recomputed with base R, independently of the Python code.
# Input: table_mr_bmi_instruments.csv written by ex11_1_bmi_covid.py (alleles already harmonised). Needs base R only.
d <- read.csv("table_mr_bmi_instruments.csv")
s <- sign(d$bx); bx <- d$bx * s; by <- d$by * s; sy <- d$sy; w <- 1 / sy^2; n <- length(bx)
ivw <- lm(by ~ 0 + bx, weights = w); b_ivw <- unname(coef(ivw)); sig <- summary(ivw)$sigma
se_fixed <- 1 / sqrt(sum(w * bx^2)); se_random <- se_fixed * max(1, sig)
egg <- lm(by ~ bx, weights = w); sm <- summary(egg); sg <- sm$sigma; se_egg <- sm$coefficients[, "Std. Error"] / sg * max(1, sg)
wald <- by / bx; se_w <- sy / abs(bx)
wmed <- function(b, wt) { o <- order(b); b <- b[o]; wt <- wt[o]; cum <- (cumsum(wt) - 0.5 * wt) / sum(wt); i <- which(cum >= 0.5)[1] - 1
  if (i < 1) return(b[1]); b[i] + (b[i + 1] - b[i]) * (0.5 - cum[i]) / (cum[i + 1] - cum[i]) }
b_wm <- wmed(wald, 1 / se_w^2)
cat(sprintf("instruments: %d\nIVW:      beta %.4f  se (fixed) %.4f  se (random) %.4f  OR %.3f\n", n, b_ivw, se_fixed, se_random, exp(b_ivw)))
cat(sprintf("MR-Egger: slope %.4f (se %.4f)  intercept %.4f (se %.4f)  OR %.3f\n", coef(egg)[2], se_egg[2], coef(egg)[1], se_egg[1], exp(coef(egg)[2])))
cat(sprintf("Weighted median: beta %.4f  OR %.3f\n", b_wm, exp(b_wm)))
Q <- sum(w * (by - b_ivw * bx)^2); cat(sprintf("Cochran Q %.1f on %d df, P = %.3f\n", Q, n - 1, pchisq(Q, n - 1, lower.tail = FALSE)))
