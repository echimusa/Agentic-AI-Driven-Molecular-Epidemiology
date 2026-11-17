# Chapter 3, Hands-on 2 (R): the per-base quality plot of the uploaded draft, drawn from the real per-base summary computed by ex3_2_read_qc.py.
# Zones as in the draft: green Q28-41 (good), orange Q20-28 (warning), red below Q20 (poor); ribbon = 25th to 75th percentile; dashed line = 10th percentile.
if (!require(ggplot2)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)
d <- read.csv("per_base_quality.csv"); d <- d[d$run %in% c("MiSeq (2014)", "NovaSeq 6000 (2024)"), ]
p <- ggplot(d, aes(x = position)) +
  annotate("rect", xmin = 0, xmax = 150, ymin = 28, ymax = 41, fill = "darkgreen", alpha = 0.15) +
  annotate("rect", xmin = 0, xmax = 150, ymin = 20, ymax = 28, fill = "orange", alpha = 0.15) +
  annotate("rect", xmin = 0, xmax = 150, ymin = 0, ymax = 20, fill = "darkred", alpha = 0.10) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = run), alpha = 0.25) +
  geom_line(aes(y = mean, color = run), linewidth = 1.1) + geom_line(aes(y = q10, color = run), linetype = "dashed", linewidth = 0.6) +
  labs(title = "Per-base sequence quality, real Salmonella runs", subtitle = "Mean, interquartile range and 10th percentile; data from the ENA", x = "Read position (bp)", y = "Phred score (Q)", color = "Run", fill = "Run") +
  scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 41)) + theme_minimal() + theme(plot.title = element_text(face = "bold", size = 13))
ggsave("figs3/fig2b_per_base_quality_R.png", p, width = 6.4, height = 4.4, dpi = 200)
cat(sprintf("%s: mean Q at base 1 = %.1f, at base 150 = %.1f\n", unique(d$run), tapply(d$mean, d$run, function(x) x[1]), tapply(d$mean, d$run, function(x) x[length(x)])), sep = "")
