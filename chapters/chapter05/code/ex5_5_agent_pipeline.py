"""
Chapter 5, Hands-on 5 (agentic) - A tool-using pipeline with checks and a human gate.
This replaces a CrewAI-style example whose agents were given no tools and whose 'autonomous' report was typed into a comment. Here each step is an explicit, logged tool call with a stated
check; a language-model agent would choose and sequence these tools, but the checks, thresholds and the final decision stay in the rules below and at the human gate. No language model or network is needed to run it.
Steps: (1) read and validate the alignment; (2) flag outlier sequences; (3) read the IQ-TREE model selection and tree; (4) report unresolved branches and weak support; (5) screen for dense clusters of substitutions; (6) test the temporal signal; (7) write a report with a
status for each check and an explicit 'needs human decision' list. Output: agent_report_<dataset>.md and agent_log_<dataset>.json.
"""
import json, re, sys, time
import numpy as np, pandas as pd
from Bio import Phylo
from phylo_utils import *
from ex5_8_recombination import screen                                          # the substitution-cluster screen of Hands-on 7
THRESH = dict(max_gap_fraction=0.05, outlier_factor=3.0, short_branch_fraction=0.5, min_support_ge80=0.7, min_r2=0.3, max_perm_p=0.01)       # thresholds are set by people and recorded in the report
def support_of(c):
    """Ultrafast bootstrap of a branch: the second number of an SH-aLRT/UFBoot label, or the label itself when only one support value is present."""
    lab = c.confidence if c.confidence is not None else c.name
    if lab is None: return None
    m = re.match(r"([\d.]+)/([\d.]+)", str(lab)); return float(m.group(2)) if m else float(lab)
def run(name, aln, prefix, dates, src=None):
    src = src or prefix
    log = []; need = []
    def step(tool, **res): log.append(dict(step=len(log) + 1, tool=tool, time=time.strftime("%H:%M:%S"), **res)); return res
    names, X = read_aln(aln); E = encode(X); L = X.shape[1]
    gaps = float((E < 0).mean()); bad = sorted(set("".join(map(str, np.unique(X)))) - set("ACGTN-RYKMSWBDHV")); step("read_alignment", taxa=len(names), sites=L, gap_or_ambiguous_fraction=round(gaps, 4), invalid_characters=bad)
    st = "pass" if gaps <= THRESH["max_gap_fraction"] and not bad else "warn"
    p, P, Q = pairwise(E); md = np.median(p, axis=1); out = [names[i] for i in np.where(md > THRESH["outlier_factor"] * np.median(md))[0]]; step("outlier_check", median_distance_to_others=round(float(np.median(md)), 5), outliers=out)
    iq = open(f"{src}.iqtree").read(); model = re.search(r"Best-fit model according to BIC: (\S+)", iq).group(1); lnl = float(re.search(r"Log-likelihood of the tree: (-?[\d.]+)", iq).group(1)); const = float(re.search(r"Number of constant sites: \d+ \(= ([\d.]+)%", iq).group(1)); step("read_iqtree", best_model_bic=model, lnL=lnl, constant_sites_pct=const)
    tr = Phylo.read(f"{src}.treefile", "newick"); internal = [c for c in tr.get_nonterminals() if c is not tr.root]; short = sum(1 for c in internal if (c.branch_length or 0) < 0.5 / L); sup = [support_of(c) for c in internal if support_of(c) is not None and (c.branch_length or 0) >= 0.5 / L]
    ge80 = float(np.mean(np.array(sup) >= 80)) if sup else float("nan"); step("branch_resolution", internal_branches=len(internal), shorter_than_one_substitution=short, share_short=round(short / len(internal), 3), resolved_support_ge80=round(ge80, 3))
    mask = screen(tr, names, E); flagged = [(names[i], int(np.where(mask[i])[0].min()) + 1, int(np.where(mask[i])[0].max()) + 1) for i in np.where(mask.any(1))[0]]; step("recombination_screen", genomes_flagged=[f[0] for f in flagged], windows=[f"{f[1]}-{f[2]}" for f in flagged])
    sig = pd.read_csv("table_root_to_tip.csv"); row = sig[(sig.dataset.str.startswith(name.split()[0])) & (sig.rooting == "best")].iloc[0]; step("temporal_signal", best_root_r2=round(float(row.r2), 3), permutation_p=float(row.permutation_p))
    checks = [("Alignment validity", st, f"gap/ambiguous {gaps:.1%}; invalid characters {bad or 'none'}"), ("Outlier sequences", "pass" if not out else "warn", f"{len(out)} flagged: {out}"),
              ("Substitution model", "info", f"BIC selects {model}; {const:.1f}% of sites are constant"), ("Tree resolution", "pass" if short / len(internal) < THRESH["short_branch_fraction"] else "warn", f"{short} of {len(internal)} internal branches are shorter than one substitution ({short / len(internal):.0%}); topology there is arbitrary"),
              ("Support", "pass" if ge80 >= THRESH["min_support_ge80"] else "warn", f"{ge80:.0%} of resolved branches have ultrafast bootstrap support of 80 or more"), ("Substitution clusters", "pass" if not flagged else "warn", "no genome with a dense cluster of private substitutions" if not flagged else f"{len(flagged)} genomes with a dense cluster of private substitutions: " + "; ".join(f"{a} (sites {b:,} to {c:,})" for a, b, c in flagged) + "; sequencing error or recombination, for a person to judge"), ("Temporal signal", "pass" if row.r2 >= THRESH["min_r2"] and row.permutation_p <= THRESH["max_perm_p"] else "warn", f"best-root R2 {row.r2:.2f}, permutation p {row.permutation_p:.3f}")]
    need = ["Confirm the thresholds above are appropriate for this question and pathogen.", "Decide whether warnings change the use of the tree (for example, do not interpret branches shorter than one substitution).", "Approve before any result is released or any alert is raised."]
    md_txt = f"# Agent report: {name}\n\nGenerated by a scripted tool pipeline; see agent_log_{prefix}.json for every call.\n\n| Check | Status | Detail |\n| --- | --- | --- |\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in checks) + "\n\n## Thresholds used (set by people)\n\n" + json.dumps(THRESH) + "\n\n## Needs a human decision\n\n" + "\n".join(f"- {x}" for x in need) + "\n"
    open(f"agent_report_{prefix}.md", "w").write(md_txt); json.dump(log, open(f"agent_log_{prefix}.json", "w"), indent=1); print(md_txt)
eb = pd.read_csv("ebola_dates.csv"); run("Ebola virus (195 genomes)", "ebola_aln.fasta", "ebola", eb, src="ebola_full"); run("H3N2 HA (60 sequences)", "h3n2_ha_aln.fasta", "h3n2", None)
