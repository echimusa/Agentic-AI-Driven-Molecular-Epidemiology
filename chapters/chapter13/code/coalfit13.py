"""Exponential-growth coalescent likelihood for a set of sampling and coalescent events, and its maximum-likelihood growth rate."""
import numpy as np
from scipy import optimize

def loglik(ev, r, logN0):
    """ev: list of (time before the last sample, +1 for a sampled lineage, -1 for a coalescence), sorted by time. N(t) = N0 exp(-r t)."""
    N0 = np.exp(logN0); k = 0; ll = 0.0; tprev = 0.0
    for t, d in ev:
        if k >= 2 and t > tprev:
            ll -= k * (k - 1) / 2 / N0 * (((np.exp(r * t) - np.exp(r * tprev)) / r) if abs(r) > 1e-9 else (t - tprev))
        if d < 0: ll += np.log(k * (k - 1) / 2) - logN0 + r * t; k -= 1
        else: k += d
        tprev = t
    return ll

def fit_r(ev, bounds=((-0.05, 0.3), (np.log(0.5), np.log(1e5))), x0=(0.05, np.log(100.0))):
    f = optimize.minimize(lambda p: -loglik(ev, p[0], p[1]), x0=list(x0), method="L-BFGS-B", bounds=list(bounds)); return f.x[0], f.x[1]

def se_r(ev, r_hat, lnN0, h=0.01):
    """Standard error of r from the curvature of the profile log-likelihood (N0 re-optimised at each r)."""
    prof = lambda r: -optimize.minimize_scalar(lambda q: -loglik(ev, r, q), bounds=(np.log(0.5), np.log(1e5)), method="bounded").fun
    d2 = (prof(r_hat + h) - 2 * prof(r_hat) + prof(r_hat - h)) / h ** 2; return 1 / np.sqrt(-d2) if d2 < 0 else np.nan
