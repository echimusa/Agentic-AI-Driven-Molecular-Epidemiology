# Chapter 7, R companion - Bayesian reconstruction with outbreaker2 on the Hagelloch measles data.
# Same public data as Exercise 1 (CRAN package `outbreaks`); timing information only.
library(outbreaks)
library(outbreaker2)

d <- measles_hagelloch_1861
onset <- as.Date(d$date_of_prodrome)
day <- as.integer(onset - min(onset))                      # days since the first prodrome

# serial interval from Exercise 1: gamma(shape 39.2, scale 0.27), discretised for lags 1..60
s <- 1:60
w <- pgamma(s + 0.5, shape = 39.2, scale = 0.27) - pgamma(s - 0.5, shape = 39.2, scale = 0.27)
w <- w / sum(w)

# onset day is treated as the infection day (delay 0), and every link is a direct one (kappa = 1)
data <- outbreaker_data(dates = day, w_dens = w, f_dens = 1)
config <- create_config(n_iter = 20000, sample_every = 20, find_import = FALSE,
                        move_kappa = FALSE, init_kappa = 1, max_kappa = 1,
                        move_pi = FALSE, init_pi = 1, move_mu = FALSE, init_mu = 1e-4,
                        move_alpha = TRUE, move_t_inf = TRUE)
set.seed(1)
res <- outbreaker(data = data, config = config)

saveRDS(res, "hagelloch_outbreaker2.rds")
keep <- res$step > 5000                                    # discard the first 5000 iterations
alpha <- res[keep, grep("^alpha_", names(res))]          # sampled infector of each case
truth <- as.integer(d$infector)                            # documented infector (NA for 4 index cases)

known <- which(!is.na(truth))
p_true <- sapply(known, function(j) mean(alpha[[j]] == truth[j], na.rm = TRUE))
modal <- sapply(known, function(j) as.integer(names(which.max(table(alpha[[j]])))))
cat(sprintf("posterior samples kept: %d\n", sum(keep)))
cat(sprintf("top-ranked infector is the documented one: %.1f%%\n", 100 * mean(modal == truth[known], na.rm = TRUE)))
cat(sprintf("mean posterior probability on the documented infector: %.1f%%\n", 100 * mean(p_true)))
cat(sprintf("cases whose documented infector has posterior > 50%%: %d of %d\n", sum(p_true > 0.5), length(known)))
