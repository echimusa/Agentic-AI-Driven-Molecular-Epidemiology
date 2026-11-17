"""
Shared scenario and tools for the cluster-triage agents of Chapter 1 (rule-table planner in ex1_6_agent_loop.py, language-model planner in ex1_7_llm_agent_loop.py).
Scenario (epi_sim.py): a known cluster of 25 sequenced cases from one simulated outbreak and nine new genomes: four from the same outbreak, two from an importation whose source lineage had evolved 14 days since the split,
two from an importation 90 days since the split, and one from the same outbreak with poor coverage. The tools are deterministic for a given case: the thinning of distances by missing sites is seeded per case, so the
result of a tool call does not depend on the order or the number of calls, which makes two agents comparable.
"""
import json
import numpy as np
from epi_sim import *
POLICY = dict(link_max_snps=2, ambiguous_max_snps=5, min_callable=0.90, max_steps=8)                  # thresholds are set by people, not by the agent
DECISIONS = ["consistent with the cluster", "ambiguous: escalate to a person", "not linked to the cluster", "not assessed: coverage below the quality gate; resequence"]
rng = np.random.default_rng(12); _ctr = [0]
_o1 = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=120, min_cases=80, mut_counter=_ctr); _idx = rng.permutation(_o1["n"])
CLUSTER = [sample_genome(_o1, i, rng.uniform(3, 10)) for i in _idx[:25]]; _same = [sample_genome(_o1, i, rng.uniform(3, 10)) for i in _idx[25:29]]
def _importation(tdiv, n):
    r = frozenset(range(10 ** 7 + _ctr[0], 10 ** 7 + _ctr[0] + rng.poisson(LAM_DEFAULT * tdiv))); _ctr[0] += 1000
    o = simulate_outbreak(rng, R=1.5, k=1.0, n_cap=40, min_cases=20, root_genome=r, mut_counter=_ctr); return [sample_genome(o, i, rng.uniform(3, 10)) for i in range(n)]
_imp14 = _importation(14, 2); _imp90 = _importation(90, 2); _poor = [sample_genome(_o1, _idx[30], 6.0)]
CASES = [("new-1", _same[0], "same outbreak", 0.99), ("new-2", _same[1], "same outbreak", 0.98), ("new-3", _imp14[0], "importation, 14 days", 0.97), ("new-4", _imp90[0], "importation, 90 days", 0.99), ("new-5", _same[2], "same outbreak", 0.97),
         ("new-6", _imp14[1], "importation, 14 days", 0.96), ("new-7", _poor[0], "same outbreak, poor coverage", 0.62), ("new-8", _imp90[1], "importation, 90 days", 0.98), ("new-9", _same[3], "same outbreak", 0.99)]
BY_NAME = {c[0]: (k, c) for k, c in enumerate(CASES)}
# ---- tools (each takes a case name and returns a dict; none changes any state)
def qc_check(name):
    c = CASES[BY_NAME[name][0]][3]; return {"callable_fraction": c, "passes_quality_gate": bool(c >= POLICY["min_callable"]), "required": POLICY["min_callable"]}
def nearest_distance(name):
    k, (_, g, _, c) = BY_NAME[name]; r = np.random.default_rng([12, k]); d = []
    for m in CLUSTER: true = len(g ^ m); d.append(int(r.binomial(true, c)) if c < 1 else true)                # missing sites hide differences
    d = np.array(d); return {"nearest_snps": int(d.min()), "cluster_members_within_2_snps": int((d <= POLICY["link_max_snps"]).sum()), "cluster_members_within_5_snps": int((d <= POLICY["ambiguous_max_snps"]).sum()), "cluster_size": len(CLUSTER)}
def decision_by_policy(nearest):
    if nearest <= POLICY["link_max_snps"]: return DECISIONS[0]
    if nearest <= POLICY["ambiguous_max_snps"]: return DECISIONS[1]
    return DECISIONS[2]
def rule_agent(name):
    """The rule-table planner of ex1_6: qc, then distance, then decide, in a fixed order; returns (decision, needs_human, steps)."""
    qc = qc_check(name)
    if not qc["passes_quality_gate"]: return DECISIONS[3], True, 2
    n = nearest_distance(name)["nearest_snps"]; d = decision_by_policy(n); return d, d == DECISIONS[1], 4
