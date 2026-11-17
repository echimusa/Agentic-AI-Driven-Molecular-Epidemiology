"""Chapter 13 helpers: generation-interval weights, Cori-style Rt, log-linear growth rate, growth rate to reproduction number."""
import numpy as np
from scipy import stats

def gamma_si(mean=15.3, sd=9.3, n=60):
    """Discretised gamma serial (or generation) interval, weights for lags 1..n days."""
    k, th = (mean / sd) ** 2, sd ** 2 / mean; edges = np.arange(0.5, n + 1.5); w = np.diff(stats.gamma.cdf(np.r_[0, edges[1:]], a=k, scale=th)); return w / w.sum()

def cori_rt(inc, w, window=14, a=1.0, b=5.0):
    """Posterior of the instantaneous reproduction number over trailing windows (Cori et al. 2013): Gamma(a + sum I, rate 1/b + sum Lambda)."""
    T = len(inc); lam = np.zeros(T)
    for t in range(1, T): s = min(t, len(w)); lam[t] = np.dot(w[:s], inc[t - 1::-1][:s])
    out = []
    for t in range(window, T):
        I, L = inc[t - window + 1:t + 1].sum(), lam[t - window + 1:t + 1].sum(); shape, rate = a + I, 1 / b + L
        out.append((t, shape / rate, *stats.gamma.ppf([0.025, 0.975], a=shape, scale=1 / rate)))
    return np.array(out)

def poisson_growth(counts, t=None):
    """Log-linear Poisson regression of counts on time: returns growth rate per day and its standard error (iteratively reweighted least squares)."""
    y = np.asarray(counts, float); t = np.arange(len(y), dtype=float) if t is None else np.asarray(t, float); X = np.column_stack([np.ones_like(t), t - t.mean()]); beta = np.array([np.log(max(y.mean(), 0.5)), 0.0])
    for _ in range(50):
        mu = np.exp(X @ beta); W = mu; z = X @ beta + (y - mu) / mu; XtW = X.T * W; new = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(new - beta)) < 1e-9: beta = new; break
        beta = new
    mu = np.exp(X @ beta); cov = np.linalg.inv((X.T * mu) @ X); phi = max(1.0, np.sum((y - mu) ** 2 / mu) / (len(y) - 2))          # quasi-Poisson scale for overdispersion
    return beta[1], np.sqrt(cov[1, 1] * phi)

def r_to_R(r, mean=15.3, sd=9.3):
    """Reproduction number from the exponential growth rate for a gamma generation interval (Wallinga and Lipsitch 2007)."""
    k, th = (mean / sd) ** 2, sd ** 2 / mean; return (1 + r * th) ** k
