"""
Chapter 6, Hands-on 4 - Time-dependent rates: why a rate measured over years should not date events over centuries (a model, not data).
Model: neutral substitutions accumulate at mu_n per genome per year; slightly deleterious mutations arise at mu_d per genome per year and are removed by purifying selection after an exponential time with rate lam per year, so a lineage of
age tau carries on average N(tau) = mu_n tau + (mu_d / lam) (1 - exp(-lam tau)) substitutions and the apparent rate is r(tau) = N(tau) / tau. The parameters are assumptions (mu_n = 10, mu_d = 40 per genome per year, lam = 0.02 per year,
a mean persistence of 50 years), chosen so that a rate measured over a decade is several times the long-term rate, the order of magnitude reported for ancient-DNA and fossil calibrations in the Adelie penguin (Subramanian et al., 2009).
Part 1 checks the closed form by simulating every mutation (arrival times, persistence times) for lineages of 1 to 1,000 years; part 2 calibrates the clock over 10 years and applies it to deeper divergences; part 3 gives the reproduction number
implied by a growth rate for gamma-distributed generation intervals (R = (1 + r theta)^a, Wallinga and Lipsitch, 2007).
Output: table_time_dependent.csv, table_time_dependent_dating.csv, table_growth_R.csv, figs6/fig_ch6_time_dependent.png.
"""
import numpy as np, pandas as pd, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
MU_N, MU_D, LAM = 10.0, 40.0, 0.02
N = lambda tau: MU_N * tau + MU_D * (1 - np.exp(-LAM * tau)) / LAM
rate = lambda tau: N(tau) / tau
rng = np.random.default_rng(3); rows = []
for tau in (1, 10, 100, 1000):
    reps = 2000; tot = np.empty(reps)
    for i in range(reps):
        k = rng.poisson(MU_D * tau); u = rng.uniform(0, tau, k); life = rng.exponential(1 / LAM, k); surv = int((u + life > tau).sum()); tot[i] = rng.poisson(MU_N * tau) + surv
    rows.append(dict(tau_years=tau, closed_form_rate=rate(tau), simulated_rate=tot.mean() / tau, simulated_sd=tot.std() / tau, ratio_to_neutral=rate(tau) / MU_N))
rows.append(dict(tau_years=10000, closed_form_rate=rate(10000), simulated_rate=np.nan, simulated_sd=np.nan, ratio_to_neutral=rate(10000) / MU_N))
T = pd.DataFrame(rows); T["relative_error_pct"] = 100 * (T.simulated_rate - T.closed_form_rate) / T.closed_form_rate; T.round(4).to_csv("table_time_dependent.csv", index=False)
r10 = rate(10); D = pd.DataFrame([dict(true_age_years=a, expected_substitutions=N(a), estimated_age_years=N(a) / r10, ratio_estimated_to_true=N(a) / r10 / a) for a in (10, 100, 1000, 10000)]); D.round(4).to_csv("table_time_dependent_dating.csv", index=False)
G = []
for td in (30, 14, 7):
    r = np.log(2) / td
    for mean in (10, 15, 20):
        a = 2.0; theta = mean / a; G.append(dict(doubling_time_days=td, growth_rate_per_day=r, mean_generation_interval_days=mean, shape=a, R=(1 + r * theta) ** a, R_exponential_interval=1 + r * mean))
pd.DataFrame(G).round(4).to_csv("table_growth_R.csv", index=False)
os.makedirs("figs6", exist_ok=True); fig, ax = plt.subplots(1, 2, figsize=(11, 4.2)); x = np.logspace(0, 4, 200)
ax[0].plot(x, rate(x) / MU_N, color="#2b6cb0", lw=2, label="closed form"); ax[0].plot(T.tau_years[:4], T.simulated_rate[:4] / MU_N, "o", color="#c05621", label="simulation"); ax[0].set(xscale="log", xlabel="Timescale of the measurement (years)", ylabel="Apparent rate / neutral rate", title="Rates depend on the timescale"); ax[0].legend(frameon=False)
a = np.logspace(1, 4, 100); ax[1].plot(a, a, color="#718096", ls=":", label="correct age"); ax[1].plot(a, N(a) / r10, color="#c05621", lw=2, label="age from a rate measured over 10 years"); ax[1].set(xscale="log", yscale="log", xlabel="True age of the divergence (years)", ylabel="Estimated age (years)", title="A short-term rate underestimates deep ages"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig("figs6/fig_ch6_time_dependent.png", dpi=180); print(T.round(3).to_string(index=False)); print(D.round(3).to_string(index=False)); print(pd.DataFrame(G).round(2).to_string(index=False))
