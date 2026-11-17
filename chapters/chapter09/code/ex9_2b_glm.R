# Chapter 9, Hands-on 2 (R companion): the adjusted odds ratios of genotypic MDR by sector, as a cross-check of the Python model.
# Input: frame_for_r.csv.gz written by ex9_2_resistance_by_sector.py. Needs base R and ggplot2.
library(ggplot2)
d <- read.csv(gzfile("frame_for_r.csv.gz"))
d$sector <- relevel(factor(d$sector), ref = "Human")
d$serovar_grp <- factor(d$serovar_grp); d$year <- factor(d$year)
fit <- glm(mdr_geno ~ sector + serovar_grp + year, data = d, family = binomial())       # model-based standard errors, as in the plain Python fit
est <- cbind(OR = exp(coef(fit)), exp(confint.default(fit)))                              # Wald intervals
out <- est[grep("^sector", rownames(est)), ]; rownames(out) <- sub("^sector", "", rownames(out)); print(round(out, 3))
res <- data.frame(sector = rownames(out), OR = out[, 1], lo = out[, 2], hi = out[, 3])
p <- ggplot(res, aes(OR, reorder(sector, OR))) + geom_vline(xintercept = 1, linetype = 2) + geom_pointrange(aes(xmin = lo, xmax = hi)) +
     scale_x_log10(breaks = c(0.5, 1, 2, 4)) + labs(x = "Adjusted odds ratio for genotypic MDR vs human (log scale)", y = NULL, title = "R check: sector and genotypic MDR") + theme_minimal(base_size = 12)
ggsave("figs9/fig5b_r_forest.png", p, width = 6.2, height = 2.8, dpi = 200)
