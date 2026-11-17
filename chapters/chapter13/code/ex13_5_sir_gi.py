"""
Chapter 13, Hands-on 5 - The SIR model of the Session 13 demonstration, and why R = exp(r T) is not enough.
A. SIR with beta = 0.3, gamma = 0.1, N = 10,000 and 10 infected at the start (R0 = beta / gamma = 3): early growth rate, doubling time, peak, final size (numerical against the final-size equation), and a control run with beta halved on day 30.
B. The same growth rate converted to R under four generation-interval shapes with the same mean (exponential, gamma with shape 2 and 4, fixed), for the SIR example and for the Ebola growth rate of Hands-on 1.
Input: r_cases.npy (Hands-on 1). Output: table_sir.csv, table_gi_shape.csv, figs13/fig5_sir.png.
"""
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False})
N, B, G, I0 = 10_000.0, 0.3, 0.1, 10.0
def sir(t, y, beta_fn): S, I, R = y; b = beta_fn(t); return [-b * S * I / N, b * S * I / N - G * I, G * I]
def run(beta_fn, T=200):
    o = solve_ivp(sir, (0, T), [N - I0, I0, 0.0], t_eval=np.arange(0, T + 0.05, 0.1), args=(beta_fn,), rtol=1e-10, atol=1e-10, method="LSODA"); return o.t, o.y        # a grid of 0.1 day, so that the peak is found to a tenth of a day
tf, (Sf, If, Rf) = run(lambda t: B)
tcf, (Scf, Icf, Rcf) = run(lambda t: B if t < 30 else B / 2)
t, S, I, R = tf[::10], Sf[::10], If[::10], Rf[::10]; tc, Sc, Ic, Rc = tcf[::10], Scf[::10], Icf[::10], Rcf[::10]
# early growth: slope of log I over the first 10 days (susceptibles are still about N)
r_est = np.polyfit(t[:11], np.log(I[:11]), 1)[0]; r_th = B - G
def final_size(R0, s0=1 - I0 / N): return brentq(lambda z: np.log(s0 / (s0 - z)) - R0 * z + 0, 1e-9, s0 - 1e-12) if False else 1 - brentq(lambda s: np.log(s0 / s) - R0 * (s0 - s), 1e-9, s0 - 1e-12)   # fraction ever infected, from ln(S0 / S_inf) = R0 (S0 - S_inf)
rows = [dict(quantity="R0 = beta / gamma", value=B / G, note="beta 0.3, gamma 0.1 per day"), dict(quantity="Early growth rate per day (log-linear slope, days 0 to 10)", value=r_est, note=f"theory beta - gamma = {r_th:.3f}"), dict(quantity="Doubling time (days) = ln 2 / r", value=np.log(2) / r_est, note=f"theory {np.log(2) / r_th:.2f}"),
        dict(quantity="Peak infected", value=If.max(), note=f"on day {tf[If.argmax()]:.1f}"), dict(quantity="Fraction ever infected, numerical", value=(N - S[-1]) / N, note="day 200"), dict(quantity="Fraction ever infected, final-size equation", value=final_size(B / G), note="ln(S0/S_inf) = R0 (S0 - S_inf)"),
        dict(quantity="Instantaneous R at the peak = R0 S / N", value=B / G * Sf[If.argmax()] / N, note="equals 1 at the peak"), dict(quantity="Control: peak infected with beta halved on day 30", value=Icf.max(), note=f"on day {tcf[Icf.argmax()]:.1f}"), dict(quantity="Control: fraction ever infected", value=(N - Sc[-1]) / N, note="day 200")]
T5 = pd.DataFrame(rows); T5.to_csv("table_sir.csv", index=False); print(T5.round(4).to_string(index=False))
# ---- B. growth rate to R under four generation-interval shapes with the same mean
def R_gamma(r, mean, shape): th = mean / shape; return (1 + r * th) ** shape
r_eb = float(np.load("r_cases.npy")[0]); sd_eb = 9.3; shape_eb = (15.3 / sd_eb) ** 2
cases = [("SIR example (r = 0.200, mean interval 10 days)", 0.2, 10.0, [1, 2, 4]), ("Ebola, Sierra Leone (r = %.4f, mean interval 15.3 days)" % r_eb, r_eb, 15.3, [1, shape_eb, 4])]
rows = []
for lab, r, m, ks in cases:
    for name, k in zip(("Exponential (shape 1)", f"Gamma, shape {ks[1]:.4g}", "Gamma, shape 4"), ks): rows.append(dict(setting=lab, interval=name, shape=k, R=R_gamma(r, m, k)))
    rows.append(dict(setting=lab, interval="Fixed (R = exp(r T))", shape=np.inf, R=np.exp(r * m)))
G5 = pd.DataFrame(rows); G5.to_csv("table_gi_shape.csv", index=False); print(G5.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.4), gridspec_kw={"width_ratios": [1.1, 0.9, 1.0]})
for y, lab, col in ((S, "Susceptible", "#2b6cb0"), (I, "Infected", "#c53030"), (R, "Recovered", "#2f855a")): ax[0].plot(t, y, color=col, lw=2, label=lab)
ax[0].plot(tc, Ic, color="#c53030", lw=1.6, ls="--", label="Infected, beta halved on day 30"); ax[0].axvline(30, color="grey", lw=0.8, ls=":"); ax[0].set(xlabel="Day", ylabel="People", title="A. SIR with R0 = 3", xlim=(0, 120)); ax[0].legend(frameon=False, fontsize=8.5)
ax[1].semilogy(t[:25], I[:25], "o", color="#c53030", ms=4); ax[1].semilogy(t[:25], I0 * np.exp(r_th * t[:25]), color="#1a202c", lw=1.4); ax[1].set(xlabel="Day", ylabel="Infected (log scale)", title="B. Early growth: slope beta - gamma"); ax[1].text(12, 12, f"r = {r_est:.3f} per day\ndoubling {np.log(2) / r_est:.2f} days", fontsize=9)
lab = ["Exponential", "Gamma\n2 or 2.7", "Gamma\n4", "Fixed"]; x = np.arange(4); a = G5[G5.setting.str.startswith("SIR")].R.values; b = G5[G5.setting.str.startswith("Ebola")].R.values
ax[2].bar(x - 0.2, a, 0.4, color="#2b6cb0", label="SIR example (r = 0.2)"); ax[2].set_ylabel("R from the growth rate (SIR example)"); ax2 = ax[2].twinx(); ax2.bar(x + 0.2, b, 0.4, color="#d69e2e", label="Ebola (r = 0.0192)"); ax2.set_ylabel("R from the growth rate (Ebola)"); ax2.set_ylim(1, 1.5); ax2.spines["top"].set_visible(False)
ax[2].set(xticks=x, xticklabels=lab, title="C. R depends on the interval's shape"); ax[2].tick_params(axis="x", labelsize=9)
h1, l1 = ax[2].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax[2].legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.5, loc="upper left")
fig.tight_layout(); fig.savefig("figs13/fig5_sir.png", dpi=200); plt.close(fig); print("figure written")
