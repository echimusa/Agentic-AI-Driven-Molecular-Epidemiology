# Chapter 10, Hands-on 2 (R companion): recompute the lead-variant P values and odds ratios from beta and SE, independently of the Python code.
# Input: qc_leads_for_r.csv written by ex10_2_qc_manhattan.py. Needs base R and ggplot2.
library(ggplot2)
d <- read.csv("qc_leads_for_r.csv")
d$z <- d$beta / d$se
d$p_wald <- 2 * pnorm(-abs(d$z))
d$log10_diff <- abs(log10(d$p_wald) - log10(d$p_reported))
d$OR <- exp(d$beta); d$lo <- exp(d$beta - 1.96 * d$se); d$hi <- exp(d$beta + 1.96 * d$se)
print(cbind(rsid = d$rsid, round(d[, c("z", "OR", "lo", "hi", "log10_diff")], 3)))
cat("largest disagreement between reported and recomputed P (log10 units):", round(max(d$log10_diff), 3), "\n")
p <- ggplot(d, aes(OR, reorder(rsid, OR))) + geom_vline(xintercept = 1, linetype = 2) + geom_pointrange(aes(xmin = lo, xmax = hi)) +
     scale_x_log10(breaks = c(0.8, 1, 1.25, 1.6)) + labs(x = "Odds ratio for hospitalised COVID-19 (95% CI, log scale)", y = NULL, title = "R check: lead variants") + theme_minimal(base_size = 12)
ggsave("figs10/fig2b_r_leads.png", p, width = 6.2, height = 3.4, dpi = 200)
