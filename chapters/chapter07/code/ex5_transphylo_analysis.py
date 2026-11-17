"""
Chapter 7, R companion 2 - score the TransPhylo runs against the simulated truth.
Reads transphylo_{runs,truth,pairs}.csv (long chains) and transphylo_short_{...}.csv (2,000-iteration chains) written by ex5_transphylo.R.
For every sampled host j, the true direct infector is a sampled host or an unsampled one (or the source). The posterior matrix gives P(i infected j directly) for sampled i.
Metrics per run (pooled over runs in the summary: links, top-1 correct, calls above 0.5 and the true ones are summed before the ratios are taken): share of hosts whose true infector was sampled, top-1 accuracy and mean posterior probability on the true infector among those hosts, precision and recall of the links with probability above 0.5,
the mean posterior probability that a host has no sampled direct infector (for hosts whose true infector was sampled, and for those whose was not), and the calibration of all pairs with probability of at least 0.01.
Writes table_transphylo_runs.csv, table_transphylo_summary.csv, table_transphylo_calibration.csv and figs3/fig12_transphylo.png.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False}); os.makedirs("figs3", exist_ok=True)
def score(prefix, label):
    R = pd.read_csv(f"{prefix}_runs.csv"); T = pd.read_csv(f"{prefix}_truth.csv"); Pp = pd.read_csv(f"{prefix}_pairs.csv"); out, cal = [], []
    for r in R.itertuples():
        ns = int(r.sampled); t = T[T.run == r.run].sort_values("host"); inf = t.infector.to_numpy()[:ns].astype(int); p = Pp[Pp.run == r.run]; M = np.zeros((ns, ns)); M[p.i.to_numpy() - 1, p.j.to_numpy() - 1] = p.prob.to_numpy()
        true_inf = np.where((inf >= 1) & (inf <= ns), inf, 0); has = true_inf > 0; links = int(has.sum())
        ok_ = [bool(M[:, j].max() > 0 and (M[:, j].argmax() + 1) == true_inf[j]) for j in range(ns) if has[j]]; top1 = np.mean(ok_) if links else np.nan; top1_correct = int(sum(ok_))
        mprob = np.mean([M[true_inf[j] - 1, j] for j in range(ns) if has[j]]) if links else np.nan
        true_mat = np.zeros((ns, ns), bool)
        for j in range(ns):
            if has[j]: true_mat[true_inf[j] - 1, j] = True
        pred = M > 0.5; prec = (pred & true_mat).sum() / pred.sum() if pred.sum() else np.nan; rec = (pred & true_mat).sum() / links if links else np.nan
        pnone = np.clip(1 - M.sum(axis=0), 0, 1)
        out.append(dict(chains=label, run=r.run, pi_true=r.pi_true, sampled=ns, hosts=r.hosts, hosts_with_sampled_infector=links, links_with_sampled_infector=links, share_sampled_infector=links / ns, top1=top1, top1_correct=top1_correct, true_called=int((pred & true_mat).sum()), pairs_called=int(pred.sum()), mean_prob_true=mprob, predicted_links=int(pred.sum()), precision_50=prec, recall_50=rec, p_no_inf_when_sampled=pnone[has].mean() if links else np.nan, p_no_inf_when_unsampled=pnone[~has].mean(), ess_log_pTTree=r.ess_log_pTTree, pi_post=r.pi_post, off_r_post=r.off_r_post, neg_post=r.neg_post, iterations=r.iterations))
        mask = M >= 0.01; cal.append(pd.DataFrame(dict(chains=label, pi_true=r.pi_true, prob=M[mask], truth=true_mat[mask])))
    return pd.DataFrame(out), pd.concat(cal)
parts = {l_: score(p_, l_) for p_, l_ in (("transphylo", "long"), ("transphylo_short", "short")) if os.path.exists(f"{p_}_runs.csv")}
bins = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0001]
for l_, (Rn_, Cal_) in parts.items():
    suf = "" if l_ == "long" else "_short"; Rn_.round(4).to_csv(f"table_transphylo{suf}_runs.csv", index=False)
    S_ = Rn_.groupby("pi_true").agg(runs=("run", "size"), sampled_mean=("sampled", "mean"), links=("links_with_sampled_infector", "sum"), share_sampled_infector=("share_sampled_infector", "mean"), top1_n=("top1_correct", "sum"), top1_runmean=("top1", "mean"), mean_prob_true=("mean_prob_true", "mean"), pairs_called=("pairs_called", "sum"), true_called=("true_called", "sum"), p_no_inf_when_sampled=("p_no_inf_when_sampled", "mean"), p_no_inf_when_unsampled=("p_no_inf_when_unsampled", "mean"), ess_median=("ess_log_pTTree", "median"), ess_min=("ess_log_pTTree", "min"), pi_post=("pi_post", "mean"), off_r_post=("off_r_post", "mean"), neg_post=("neg_post", "mean")).reset_index()
    S_["pooled_top1"] = S_.top1_n / S_.links; S_["precision"] = S_.true_called / S_.pairs_called; S_["recall"] = S_.true_called / S_.links; S_["top1"] = S_.pooled_top1; S_["precision_50"] = S_.precision; S_["recall_50"] = S_.recall; S_["p_none_if_sampled"] = S_.p_no_inf_when_sampled; S_["p_none_if_unsampled"] = S_.p_no_inf_when_unsampled
    S_.round(4).to_csv(f"table_transphylo{suf}_summary.csv", index=False)
    Cal_["bin"] = pd.cut(Cal_.prob, bins, right=False); C_ = Cal_.groupby(["pi_true", "bin"], observed=True).agg(n=("prob", "size"), mean_prob=("prob", "mean"), observed=("truth", "mean")).reset_index(); C_["bin"] = C_.bin.astype(str); C_.round(4).to_csv(f"table_transphylo{suf}_calibration.csv", index=False)
    print(l_, "chains"); print(S_.round(2).to_string(index=False))
Rn = pd.concat([p[0] for p in parts.values()]); Cal = pd.concat([p[1] for p in parts.values()]); C = pd.read_csv("table_transphylo_calibration.csv") if "long" in parts else pd.read_csv("table_transphylo_short_calibration.csv"); C["chains"] = "long" if "long" in parts else "short"
# ------------------------------------------------------------- figure (long chains)
lab = "long" if "long" in parts else "short"; L = parts[lab][0]
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2))
for k, (col, name) in enumerate([("share_sampled_infector", "True infector was sampled"), ("recall_50", "Links found (prob > 0.5)"), ("precision_50", "Links correct (prob > 0.5)")]):
    for pi, c in ((0.25, "#e4572e"), (0.75, "#1b998b")):
        v = L[L.pi_true == pi][col].dropna(); ax[0].scatter([k + (0.12 if pi > 0.5 else -0.12)] * len(v), 100 * v, color=c, s=30, alpha=0.8, label=f"pi = {pi}" if k == 0 else None)
ax[0].set(xticks=[0, 1, 2], xticklabels=["True infector\nsampled", "Links found", "Links correct"], ylabel="Per cent", ylim=(-3, 103), title="Sampling fraction and accuracy"); ax[0].legend(frameon=False)
c_ = parts[lab][1].assign(bin=lambda x: pd.cut(x.prob, bins, right=False)).groupby("bin", observed=True).agg(mean_prob=("prob", "mean"), observed=("truth", "mean")).reset_index(); ax[1].plot([0, 1], [0, 1], color="#718096", ls=":"); ax[1].plot(c_.mean_prob, c_.observed, "o-", color="#2b6cb0"); ax[1].set(xlabel="Posterior probability of a direct link", ylabel="Share of such pairs that are true", title="Calibration of the pairwise probabilities")
r0 = L.iloc[0]; T0 = pd.read_csv("transphylo.csv".replace("transphylo.csv", "transphylo_truth.csv" if lab == "long" else "transphylo_short_truth.csv")); P0 = pd.read_csv("transphylo_pairs.csv" if lab == "long" else "transphylo_short_pairs.csv")
ns = int(r0.sampled); p = P0[P0.run == r0.run]; M = np.zeros((ns, ns)); M[p.i - 1, p.j - 1] = p.prob; im = ax[2].imshow(M, cmap="Blues", vmin=0, vmax=1, origin="lower"); t = T0[T0.run == r0.run].sort_values("host").infector.to_numpy()[:ns].astype(int)
for j in range(ns):
    if 1 <= t[j] <= ns: ax[2].plot(j, t[j] - 1, "o", mfc="none", mec="#e4572e", ms=9)
ax[2].set(xlabel="Infectee j", ylabel="Infector i", title="Posterior P(i infected j) and the true links"); fig.colorbar(im, ax=ax[2], fraction=0.046); fig.tight_layout(); fig.savefig("figs3/fig12_transphylo.png", dpi=180)
