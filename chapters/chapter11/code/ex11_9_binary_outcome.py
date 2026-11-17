"""
Chapter 11, Hands-on 9 - a binary outcome: testing a causal effect is easier than estimating one (Didelez and Sheehan 2007).
A: non-collapsibility. X = 0.5 G + U + noise with G a standardised variant and U an unmeasured confounder; Y is binary with logit P(Y = 1 | X, U) = a + b X + c U, b = 0.5, c = 1. The conditional odds ratio exp(b) = 1.65 is a property of the model
   given U; the causal odds ratio for raising X by one unit is a ratio of population risks, P(Y=1 | do(X=x+1)) / P(Y=0 | do(X=x+1)) over the same at x, which averages over U and is smaller. The script computes it by numerical integration.
   It then estimates by the ratio of the logistic G-Y coefficient to the linear G-X coefficient (the Wald ratio on the log-odds scale) in 500 replicates of 100,000 people, and by ordinary logistic regression of Y on X.
B: testing. With b = 0 (X has no effect) and the same confounding, the share of replicates that reject at 5% for (i) logistic regression of Y on X and (ii) the test that G is associated with Y; the second keeps its size.
C: bounds. With G, X and Y all binary, the core instrument conditions alone give bounds on the average causal effect (Balke and Pearl 1997; Didelez and Sheehan 2007, Section 6.2). The script computes the eight lower and eight upper expressions
   for a simulated population and checks the bounds against a linear programme over the 16 response types, and against the true effect.
Output: Figure 8, table_binary_estimation.csv, table_binary_testing.csv, table_bounds.csv.   Needs numpy, pandas, scipy, matplotlib only.
"""
import os, itertools, numpy as np, pandas as pd
from scipy import stats, optimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs11", exist_ok=True)
rng = np.random.default_rng(13); N, REPS = 100_000, 500; B, C, A0, GAM = 0.5, 1.0, -2.0, 0.5
expit = lambda z: 1 / (1 + np.exp(-z))
def fit_logit(y, X, iters=30):
    X1 = np.column_stack([np.ones(len(y)), X]); beta = np.zeros(X1.shape[1])
    for _ in range(iters):
        p = expit(X1 @ beta); W = p * (1 - p); H = X1.T @ (X1 * W[:, None]); step = np.linalg.solve(H, X1.T @ (y - p)); beta += step
        if np.max(np.abs(step)) < 1e-9: break
    se = np.sqrt(np.diag(np.linalg.inv(H))); return beta[1:], se[1:]
# ---------------------------------------------------------------- A: the causal odds ratio by integration over U
u, wu = np.polynomial.hermite_e.hermegauss(80); wu = wu / wu.sum()           # N(0,1) quadrature
def risk_do(x, a, b): return float(np.sum(wu * expit(a + b * x + C * u)))
def causal_or(b, a=A0, x0=0.0, dx=1.0): r1, r0 = risk_do(x0 + dx, a, b), risk_do(x0, a, b); return (r1 / (1 - r1)) / (r0 / (1 - r0))
true_cor = causal_or(B); print(f"conditional OR per unit of X: {np.exp(B):.3f}; causal OR per unit of X: {true_cor:.3f}; causal log OR {np.log(true_cor):.3f}")
def simulate(b):
    G = rng.binomial(2, 0.3, N).astype(float); G = (G - 0.6) / np.sqrt(2 * 0.3 * 0.7); U = rng.normal(size=N); X = GAM * G + U + rng.normal(size=N); Y = (rng.random(N) < expit(A0 + b * X + C * U)).astype(float); return G, X, Y
est = {"obs": [], "iv": [], "iv_se": []}; rej = {"obs": [], "gy": []}
for _ in range(REPS):
    G, X, Y = simulate(B); o, _ = fit_logit(Y, X); est["obs"].append(o[0]); gy, gy_se = fit_logit(Y, G); gx = np.cov(G, X)[0, 1] / G.var(); est["iv"].append(gy[0] / gx); est["iv_se"].append(gy_se[0] / abs(gx))
iv, ob = np.array(est["iv"]), np.array(est["obs"]); E = pd.DataFrame([dict(quantity="Conditional log odds ratio (model parameter b)", value=B), dict(quantity="Causal log odds ratio per unit of X (integrating over U)", value=np.log(true_cor)),
     dict(quantity="Observational logistic regression, mean", value=ob.mean()), dict(quantity="IV estimate (logistic G-Y over linear G-X), mean", value=iv.mean()), dict(quantity="IV estimate, standard deviation over replicates", value=iv.std()), dict(quantity="IV estimate, mean model standard error", value=np.mean(est["iv_se"]))])
E.to_csv("table_binary_estimation.csv", index=False); print(E.round(4).to_string(index=False))
for _ in range(REPS):
    G, X, Y = simulate(0.0); o, os_ = fit_logit(Y, X); rej["obs"].append(abs(o[0] / os_[0]) > 1.96); g, gs = fit_logit(Y, G); rej["gy"].append(abs(g[0] / gs[0]) > 1.96)
T = pd.DataFrame([dict(test="Logistic regression of Y on X (confounded)", rejects_at_5pct=float(np.mean(rej["obs"])), true_effect=0.0), dict(test="Test that G is associated with Y (instrument test)", rejects_at_5pct=float(np.mean(rej["gy"])), true_effect=0.0)]); T.to_csv("table_binary_testing.csv", index=False); print(T.round(3).to_string(index=False))
# ---------------------------------------------------------------- C: Balke-Pearl bounds, and a linear programme over the 16 response types
def bounds(p):                                   # p[(y, x, z)] = P(Y = y, X = x | Z = z)
    P = lambda y, x, z: p[(y, x, z)]
    lo = [P(1, 1, 1) + P(0, 0, 0) - 1, P(1, 1, 0) + P(0, 0, 1) - 1, P(1, 1, 0) - P(1, 1, 1) - P(1, 0, 1) - P(0, 1, 0) - P(1, 0, 0), P(1, 1, 1) - P(1, 1, 0) - P(1, 0, 0) - P(0, 1, 1) - P(1, 0, 1),
          -P(0, 1, 1) - P(1, 0, 1), -P(0, 1, 0) - P(1, 0, 0), P(0, 0, 1) - P(0, 1, 1) - P(1, 0, 1) - P(0, 1, 0) - P(0, 0, 0), P(0, 0, 0) - P(0, 1, 0) - P(1, 0, 0) - P(0, 1, 1) - P(0, 0, 1)]
    hi = [1 - P(0, 1, 1) - P(1, 0, 0), 1 - P(0, 1, 0) - P(1, 0, 1), -P(0, 1, 0) + P(0, 1, 1) + P(0, 0, 1) + P(1, 1, 0) + P(0, 0, 0), -P(0, 1, 1) + P(0, 1, 0) + P(0, 0, 0) + P(1, 1, 1) + P(0, 0, 1),
          P(1, 1, 1) + P(0, 0, 1), P(1, 1, 0) + P(0, 0, 0), -P(1, 0, 1) + P(1, 1, 1) + P(0, 0, 1) + P(1, 1, 0) + P(1, 0, 0), -P(1, 0, 0) + P(1, 1, 0) + P(0, 0, 0) + P(1, 1, 1) + P(1, 0, 1)]
    return max(lo), min(hi)
types = list(itertools.product([0, 1], repeat=2) + [()] if False else itertools.product(range(2), range(2), range(2), range(2)))      # (x under z=0, x under z=1, y under x=0, y under x=1)
def lp_bounds(p):
    A, b = [], []
    for z in (0, 1):
        for x in (0, 1):
            for y in (0, 1):
                A.append([1.0 if (t[z] == x and t[2 + x] == y) else 0.0 for t in types]); b.append(p[(y, x, z)])
    ace = np.array([t[3] - t[2] for t in types], float); out = []
    for sgn in (1, -1): r = optimize.linprog(sgn * ace, A_eq=np.array(A), b_eq=np.array(b), bounds=(0, None), method="highs"); out.append(sgn * r.fun)
    return out[0], out[1]
def make_population(q):                          # q: probabilities of the 16 response types
    p = {}; ace = float(sum(qi * (t[3] - t[2]) for qi, t in zip(q, types)))
    for z in (0, 1):
        for x in (0, 1):
            for y in (0, 1): p[(y, x, z)] = float(sum(qi for qi, t in zip(q, types) if t[z] == x and t[2 + x] == y))
    return p, ace
rows = []
for k in range(8):
    q = rng.dirichlet(np.ones(16) * (0.5 if k < 4 else 3.0)); p, ace = make_population(q); lo, hi = bounds(p); llo, lhi = lp_bounds(p)
    rows.append(dict(population=k + 1, true_ace=ace, formula_lower=lo, formula_upper=hi, lp_lower=llo, lp_upper=lhi, width=hi - lo, inside=bool(lo - 1e-9 <= ace <= hi + 1e-9)))
Bd = pd.DataFrame(rows); Bd.to_csv("table_bounds.csv", index=False); print(Bd.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.2)); xs = np.linspace(-1.2, 1.2, 25)
ax[0].plot(xs, [np.log(causal_or(b)) for b in xs], color="#E4572E", lw=2, label="causal log odds ratio"); ax[0].plot(xs, xs, color="#a0aec0", lw=1.4, ls="--", label="conditional (model) log odds ratio"); ax[0].set(xlabel="Model parameter b", ylabel="Log odds ratio per unit of X", title="The two odds ratios differ"); ax[0].legend(frameon=False, fontsize=8.5, loc="upper left")
ax[1].hist(iv, bins=24, color="#3D7EA6", alpha=0.85, label="IV estimates"); ax[1].axvline(np.log(true_cor), color="#E4572E", lw=2, label="causal log OR"); ax[1].axvline(B, color="#a0aec0", lw=1.6, ls="--", label="conditional log OR"); ax[1].axvline(ob.mean(), color="#2E294E", lw=1.6, label="observational"); ax[1].set(xlabel="Estimate (log odds ratio)", ylabel="Replicates", title="Where each estimate lands"); ax[1].legend(frameon=False, fontsize=8, loc="upper left")
for i, r in Bd.iterrows(): ax[2].plot([r.formula_lower, r.formula_upper], [i, i], color="#3D7EA6", lw=5, solid_capstyle="butt"); ax[2].plot(r.true_ace, i, "o", color="#E4572E")
ax[2].axvline(0, color="#a0aec0", lw=0.8); ax[2].set(xlabel="Average causal effect on risk (orange: truth)", ylabel="Simulated population", title="The core assumptions alone give wide bounds"); ax[2].set_yticks(range(len(Bd))); ax[2].set_yticklabels(Bd.population)
fig.tight_layout(); fig.savefig("figs11/fig8_binary_outcome.png", dpi=200); print("figure written")
