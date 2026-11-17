"""
Chapter 2, Hands-on 1 - Nested case-control sampling from a cohort: bias, efficiency and cost.
Starts from the uploaded Chapter_2_1.py (cohort of 10,000, a biomarker with true log hazard 0.4 per SD, 1:2 matched controls) and
  (a) writes the control selection in vectorised form,
  (b) contrasts the uploaded rule (controls must never become cases) with risk-set (incidence density) sampling, where anyone still at risk at the case's diagnosis may be a control,
  (c) compares both with the full-cohort Cox estimate across 1, 2, 4 and 8 controls per case over 200 simulated cohorts.
Output: table_ncc_design.csv, table_infert_clogit.csv, Figure 1.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.discrete.conditional_models import ConditionalLogit
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
N, BETA, REPS, MS = 10000, 0.4, 200, [1, 2, 4, 8]

def make_cohort(rng):
    x = rng.normal(0, 1, N); p = 1 / (1 + np.exp(-(-3.5 + BETA * x))); case = rng.uniform(0, 1, N) < p
    day = np.where(case, rng.integers(1, 365, N), 365); return x, case.astype(int), day            # day = diagnosis day, or 365 for people still disease-free at the end

def select(x, case, day, M, rng, risk_set):
    """One matched set per case. risk_set=True: controls are anyone not yet diagnosed at the case's day (they may become cases later).
    risk_set=False: the uploaded rule, controls must be disease-free for the whole follow-up."""
    ids = np.where(case == 1)[0]; rows = []; sets = 0
    for i in ids:
        t = day[i]; ok = (day > t) if risk_set else ((case == 0) & (day >= t)); ok[i] = False; pool = np.flatnonzero(ok)
        if len(pool) < M: continue
        ctrl = rng.choice(pool, M, replace=False); sets += 1; rows += [(sets, x[i], 1)] + [(sets, x[c], 0) for c in ctrl]
    return np.array(rows)

def clogit_beta(arr):
    r = ConditionalLogit(arr[:, 2], arr[:, 1:2], groups=arr[:, 0]).fit(disp=0); return float(r.params[0]), float(r.bse[0])

rng = np.random.default_rng(42); out = []
for rep in range(REPS):
    x, case, day = make_cohort(rng); cox = PHReg(day.astype(float), x[:, None], status=case).fit(disp=0); b_c = float(cox.params[0]); se_c = float(cox.bse[0])
    out.append(dict(rep=rep, design="Full cohort (Cox)", M=0, beta=b_c, se=se_c, measured=N, cases=int(case.sum())))
    for M in MS:
        for rs, name in [(True, "Risk-set sampling"), (False, "Uploaded rule (never-cases)")]:
            arr = select(x, case, day, M, rng, rs); b, se = clogit_beta(arr); out.append(dict(rep=rep, design=name, M=M, beta=b, se=se, measured=len(arr), cases=int(case.sum())))
S = pd.DataFrame(out)
base = S[S.design == "Full cohort (Cox)"]; var_cohort = base.beta.var()
T = S.groupby(["design", "M"]).agg(mean_beta=("beta", "mean"), sd_beta=("beta", "std"), mean_se=("se", "mean"), measured=("measured", "mean"), cases=("cases", "mean")).reset_index()
T["bias"] = T.mean_beta - BETA; T["odds_ratio"] = np.exp(T.mean_beta); T["efficiency_vs_cohort"] = var_cohort / T.sd_beta ** 2; T["share_measured_pct"] = 100 * T.measured / N; T["coverage95"] = [((S[(S.design == d) & (S.M == m)].beta - BETA).abs() < 1.96 * S[(S.design == d) & (S.M == m)].se).mean() * 100 for d, m in zip(T.design, T.M)]
T.round(4).to_csv("table_ncc_design.csv", index=False); print(f"true log hazard per SD = {BETA} (HR {np.exp(BETA):.3f}); about {S.cases.mean():.0f} cases per cohort; {REPS} simulated cohorts"); print(T.round(3).to_string(index=False))
# a public matched case-control data set: the 'infert' data of base R (Trichopoulos et al. 1976; 83 cases, 1-2 matched controls each)
inf = pd.read_csv("infert.csv"); rows = []
for name, cols in [("spontaneous", ["spontaneous"]), ("induced", ["induced"]), ("both", ["spontaneous", "induced"])]:
    r = ConditionalLogit(inf.case.values, inf[cols].values.astype(float), groups=inf.stratum.values).fit(disp=0)
    for c, b, se in zip(cols, r.params, r.bse): rows.append(dict(model=name, term=c, log_or=b, se=se, odds_ratio=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se)))
I = pd.DataFrame(rows); I.round(4).to_csv("table_infert_clogit.csv", index=False); print(f"infert: {len(inf)} women, {inf.case.sum()} cases, {inf.stratum.nunique()} matched sets"); print(I.round(3).to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.6), gridspec_kw={"width_ratios": [1.2, 1, 1]})
labs = ["Cohort"] + [f"1:{m}" for m in MS]; data_rs = [base.beta] + [S[(S.design == "Risk-set sampling") & (S.M == m)].beta for m in MS]; data_up = [base.beta] + [S[(S.design == "Uploaded rule (never-cases)") & (S.M == m)].beta for m in MS]
pos = np.arange(len(labs)); b1 = ax[0].boxplot(data_rs, positions=pos - 0.18, widths=0.3, patch_artist=True, showfliers=False); b2 = ax[0].boxplot(data_up, positions=pos + 0.18, widths=0.3, patch_artist=True, showfliers=False)
for b in b1["boxes"]: b.set(facecolor="#2b6cb0", alpha=0.6)
for b in b2["boxes"]: b.set(facecolor="#dd6b20", alpha=0.6)
ax[0].axhline(BETA, color="black", ls="--", lw=1); ax[0].set_xticks(pos); ax[0].set_xticklabels(labs); ax[0].set(xlabel="Controls per case", ylabel="Estimated log hazard ratio per SD", title="Estimates over 200 cohorts (true value dashed)")
ax[0].legend([b1["boxes"][1], b2["boxes"][1]], ["Risk-set sampling", "Uploaded rule"], frameon=False, fontsize=8.5, loc="upper right")
for name, c in [("Risk-set sampling", "#2b6cb0"), ("Uploaded rule (never-cases)", "#dd6b20")]: x_ = T[T.design == name]; ax[1].plot(x_.M, x_.efficiency_vs_cohort, "o-", color=c, label=name)
mm = np.arange(1, 9); ax[1].plot(mm, mm / (mm + 1), "k:", label="M / (M + 1)"); ax[1].set(xlabel="Controls per case (M)", ylabel="Efficiency relative to the full cohort", title="What each extra control buys", ylim=(0, 1.15)); ax[1].legend(frameon=False, fontsize=8.5, loc="lower right")
x_ = T[T.design == "Risk-set sampling"]; ax[2].plot(x_.share_measured_pct, x_.efficiency_vs_cohort, "o-", color="#2b6cb0")
for r in x_.itertuples(): ax[2].annotate(f"1:{int(r.M)}", (r.share_measured_pct, r.efficiency_vs_cohort), xytext=(5, -12), textcoords="offset points", fontsize=8.5)
ax[2].set(xlabel="Share of the cohort measured (%)", ylabel="Efficiency relative to the full cohort", title="Efficiency against assay cost", ylim=(0, 1.15)); fig.tight_layout(); fig.savefig("figs2/fig1_ncc_design.png", dpi=200); plt.close(fig); print("figure written")
