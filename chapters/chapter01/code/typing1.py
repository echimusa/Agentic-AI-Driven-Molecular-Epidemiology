"""Chapter 1 helpers: Simpson's index of diversity (Hunter and Gaston 1988), Wallace and adjusted Wallace coefficients (Carrico 2006, Severiano 2011), and confidence intervals by subsampling (boot_ci, the bootstrap with replacement, is kept for comparison only: it duplicates isolates, and a duplicate always shares every type, so it inflates Wallace coefficients)."""
import numpy as np, pandas as pd

def simpson_D(labels):
    """Hunter-Gaston discriminatory index: probability that two isolates drawn at random have different types."""
    n = pd.Series(labels).value_counts().values.astype(float); N = n.sum(); return 1 - (n * (n - 1)).sum() / (N * (N - 1))

def _pairs(series):
    return (series * (series - 1) / 2).sum()

def wallace(a, b):
    """W(A->B): probability that two isolates sharing type A also share type B. Also returns the adjusted Wallace (chance-corrected)."""
    a, b = pd.Series(np.asarray(a)), pd.Series(np.asarray(b)); both = pd.DataFrame({"a": a, "b": b}).groupby(["a", "b"]).size().values.astype(float)
    pa = _pairs(a.value_counts().values.astype(float)); w = _pairs(both) / pa if pa > 0 else np.nan
    sb = 1 - simpson_D(b); aw = (w - sb) / (1 - sb) if sb < 1 else np.nan; return w, aw

def boot_ci(df, fn, reps=300, seed=1):
    rng = np.random.default_rng(seed); n = len(df); vals = [fn(df.iloc[rng.integers(0, n, n)]) for _ in range(reps)]; return np.nanpercentile(np.array(vals), [2.5, 97.5], axis=0)

def subsample_ci(df, fn, reps=400, frac=0.5, seed=1, z=1.96):
    """Delete-half subsampling interval: draw half of the isolates WITHOUT replacement, compute the statistic, and scale the spread of those half-sample estimates by sqrt(m / (N - m)) to the full sample size.
    The interval is estimate +/- z * scaled SD, so it is centred on the full-sample estimate and is not distorted by duplicated isolates."""
    rng = np.random.default_rng(seed); n = len(df); m = int(n * frac); est = fn(df)
    vals = np.array([fn(df.iloc[rng.choice(n, m, replace=False)]) for _ in range(reps)], dtype=float); se = np.sqrt(m / (n - m)) * np.nanstd(vals, axis=0, ddof=1)
    return est - z * se, est + z * se
