# Chapter 2, Hands-on 1 (R check): a conditional logistic regression written from the likelihood in base R, on the matched case-control data 'infert'.
# Each matched set has one case and one or more controls; the set's log-likelihood is  x_case * b - log( sum over the set of exp(x * b) ).
d <- datasets::infert
clog <- function(vars) {
  X <- as.matrix(d[, vars, drop = FALSE]); sets <- split(seq_len(nrow(d)), d$stratum)
  nll <- function(b) -sum(sapply(sets, function(ix) { eta <- drop(X[ix, , drop = FALSE] %*% b); sum(eta[d$case[ix] == 1]) - log(sum(exp(eta))) }))
  fit <- optim(rep(0, length(vars)), nll, method = "BFGS", hessian = TRUE); se <- sqrt(diag(solve(fit$hessian)))
  data.frame(term = vars, log_or = fit$par, se = se, odds_ratio = exp(fit$par), lo = exp(fit$par - 1.96 * se), hi = exp(fit$par + 1.96 * se))
}
cat(sprintf("%d women; %d cases; %d matched sets\n", nrow(d), sum(d$case), length(unique(d$stratum))))
print(rbind(cbind(model = "spontaneous", clog("spontaneous")), cbind(model = "induced", clog("induced")), cbind(model = "both", clog(c("spontaneous", "induced")))), digits = 4, row.names = FALSE)
