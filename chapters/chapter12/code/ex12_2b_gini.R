# Chapter 12, Hands-on 2 (R companion): the Lorenz curve and Gini coefficient of isolates over population, recomputed in base R.
# Input: table_country_equity.csv written by ex12_2_equity.py. Needs base R only.
d <- read.csv("table_country_equity.csv")
o <- order(d$isolates / d$pop); iso <- as.numeric(d$isolates[o]); pop <- as.numeric(d$pop[o])
L <- c(0, cumsum(iso) / sum(iso)); P <- c(0, cumsum(pop) / sum(pop))
gini <- 1 - sum((P[-1] - P[-length(P)]) * (L[-1] + L[-length(L)]))
cat(sprintf("economies: %d | economies with at least one isolate: %d | Gini: %.3f\n", nrow(d), sum(d$isolates > 0), gini))
top3 <- sort(d$isolates, decreasing = TRUE)[1:3]; cat(sprintf("top three economies hold %.1f%% of isolates\n", 100 * sum(top3) / sum(d$isolates)))
png("figs12/fig2b_r_lorenz.png", width = 1200, height = 1000, res = 200); plot(P, L, type = "l", lwd = 2, col = "#2b6cb0", xlab = "Cumulative share of population", ylab = "Cumulative share of isolates", main = "R check: Lorenz curve"); abline(0, 1, lty = 2); dev.off()
