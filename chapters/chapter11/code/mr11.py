"""
Chapter 11 helpers: allele harmonisation and the core two-sample Mendelian randomisation estimators, written with numpy only.
Inputs are per-variant effects on the exposure (bx, sx) and on the outcome (by, sy), aligned to the same effect allele.
"""
import numpy as np
from scipy import stats

def harmonise(e_a1, e_a2, e_eaf, o_a1, o_a2, o_eaf, palindrome_window=(0.42, 0.58)):
    """Return (keep mask, sign) so that outcome effects can be aligned to the exposure effect allele. Palindromic variants are
    aligned by allele frequency and dropped when the frequency is too close to 0.5 to tell the strands apart."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    e_a1, e_a2, o_a1, o_a2 = (np.char.upper(np.asarray([str(v) for v in x], dtype="U12")) for x in (e_a1, e_a2, o_a1, o_a2))
    same = (e_a1 == o_a1) & (e_a2 == o_a2); flip = (e_a1 == o_a2) & (e_a2 == o_a1)
    c = lambda x: np.array([comp.get(v, "?") for v in x]); comp_same = (e_a1 == c(o_a1)) & (e_a2 == c(o_a2)); comp_flip = (e_a1 == c(o_a2)) & (e_a2 == c(o_a1))
    pal = (e_a1 == c(e_a2)); ok = same | flip | comp_same | comp_flip; sign = np.where(same | comp_same, 1, -1)
    amb = pal & (((e_eaf > palindrome_window[0]) & (e_eaf < palindrome_window[1])) | ((o_eaf > palindrome_window[0]) & (o_eaf < palindrome_window[1])) | np.isnan(e_eaf) | np.isnan(o_eaf))
    # palindromic and not ambiguous: align by frequency (same side of 0.5 means same orientation)
    pal_sign = np.where((e_eaf > 0.5) == (o_eaf > 0.5), 1, -1); sign = np.where(pal, pal_sign, sign)
    return ok & ~amb, sign

def wald(bx, by, sy): return by / bx, sy / np.abs(bx)

def ivw(bx, by, sy):
    w = bx ** 2 / sy ** 2; b = np.sum(bx * by / sy ** 2) / np.sum(w); se_f = 1 / np.sqrt(np.sum(w))
    Q = np.sum((by - b * bx) ** 2 / sy ** 2); df = len(bx) - 1; phi = max(1.0, Q / df) if df > 0 else 1.0
    return dict(beta=b, se_fixed=se_f, se_random=se_f * np.sqrt(phi), Q=Q, df=df, Q_p=stats.chi2.sf(Q, df) if df > 0 else np.nan, I2=max(0.0, (Q - df) / Q) if Q > 0 else 0.0)

def egger(bx, by, sy):
    """MR-Egger: weighted regression of outcome effects on exposure effects with a free intercept, after orienting variants so that bx > 0."""
    s = np.sign(bx); bx_, by_ = bx * s, by * s; w = 1 / sy ** 2; X = np.column_stack([np.ones_like(bx_), bx_]); W = np.diag(w) if len(bx_) < 2000 else None
    XtW = X.T * w; cov = np.linalg.inv(XtW @ X); beta = cov @ (XtW @ by_); res = by_ - X @ beta; df = len(bx_) - 2; sigma2 = np.sum(w * res ** 2) / df; sig = max(1.0, np.sqrt(sigma2))
    se = np.sqrt(np.diag(cov)) * sig; t = beta / se; p = 2 * stats.t.sf(np.abs(t), df)
    return dict(intercept=beta[0], se_intercept=se[0], p_intercept=p[0], slope=beta[1], se_slope=se[1], p_slope=p[1], sigma=np.sqrt(sigma2))

def weighted_median(b, se, nboot=1000, seed=1):
    """Weighted median of the Wald ratios (Bowden 2016), with a parametric bootstrap standard error."""
    def wm(bb, ww):
        o = np.argsort(bb); bb, ww = bb[o], ww[o]; cum = np.cumsum(ww) - 0.5 * ww; cum /= ww.sum(); i = np.searchsorted(cum, 0.5) - 1
        if i < 0: return bb[0]
        if i >= len(bb) - 1: return bb[-1]
        return bb[i] + (bb[i + 1] - bb[i]) * (0.5 - cum[i]) / (cum[i + 1] - cum[i])
    w = 1 / se ** 2; est = wm(b, w); rng = np.random.default_rng(seed); bs = [wm(rng.normal(b, se), w) for _ in range(nboot)]
    return dict(beta=est, se=np.std(bs, ddof=1))

def leave_one_out(bx, by, sy): return np.array([ivw(np.delete(bx, i), np.delete(by, i), np.delete(sy, i))["beta"] for i in range(len(bx))])

def fstat(bx, sx): return (bx / sx) ** 2

def i2gx(bx, sx):
    w = 1 / sx ** 2; m = np.sum(w * bx) / np.sum(w); Q = np.sum(w * (bx - m) ** 2); df = len(bx) - 1; return max(0.0, (Q - df) / Q)

# ---------------------------------------------------------------- added for Hands-on 6 to 10
def weighted_mode(b, se, phi=1.0, grid=4000, nboot=0, seed=1):
    """Weighted mode of the ratio estimates (Hartwig, Davey Smith and Bowden 2017): the peak of a kernel density of the Wald ratios, each weighted by its inverse variance.
    Bandwidth 0.9 min(sd, mad) n^(-1/5) times phi, as in the MendelianRandomization package. nboot > 0 adds a parametric bootstrap standard error."""
    b = np.asarray(b, float); se = np.asarray(se, float)
    def mode(bb, ss):
        w = 1 / ss ** 2; w = w / w.sum(); mad = 1.4826 * np.median(np.abs(bb - np.median(bb))); h = phi * 0.9 * min(np.std(bb, ddof=1), mad) / len(bb) ** 0.2
        if h <= 0: return float(np.median(bb))
        xs = np.linspace(bb.min() - 3 * h, bb.max() + 3 * h, grid); d = (w[None, :] * np.exp(-0.5 * ((xs[:, None] - bb[None, :]) / h) ** 2)).sum(axis=1)
        return float(xs[np.argmax(d)])
    est = mode(b, se); out = dict(beta=est)
    if nboot:
        rng = np.random.default_rng(seed); out["se"] = float(np.std([mode(rng.normal(b, se), se) for _ in range(nboot)], ddof=1))
    return out

def mvmr_ivw(BX, by, sy):
    """Multivariable IVW (Sanderson et al. 2019): weighted regression of the outcome effects on the effects on every exposure, through the origin. BX is variants x exposures."""
    w = 1 / sy ** 2; XtW = BX.T * w; cov = np.linalg.inv(XtW @ BX); beta = cov @ (XtW @ by); res = by - BX @ beta; df = len(by) - BX.shape[1]; Q = np.sum(w * res ** 2); phi = max(1.0, Q / df)
    return dict(beta=beta, se=np.sqrt(np.diag(cov) * phi), Q=Q, df=df)

def r2_from_z(z, n): return z ** 2 / (z ** 2 + n - 2)          # variance explained by a variant, from its Z statistic and the sample size (Steiger)
