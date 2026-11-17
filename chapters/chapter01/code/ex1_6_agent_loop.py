"""
Chapter 1, agentic AI foundations - a minimal agent loop for outbreak-cluster triage, with a human gate (rule-table planner).
The loop has the parts of any agent: a goal, a planner that chooses the next tool from the current state, tools that act, observations that update the state, checks that can stop the loop, a step budget, and a rule that
sends uncertain cases to a person. The planner here is a small rule table so the example runs without a language model or network; ex1_7_llm_agent_loop.py replaces the planner with a language model and keeps the tools,
the policy and the guardrails. Scenario and tools: triage_core.py (a known cluster of 25 sequenced cases; nine new genomes).
Output: table_agent_trace.csv, agent_trace.txt.
"""
import json
import numpy as np, pandas as pd
from triage_core import *
def draft_note(state, decision): return f"{state['name']}: {decision}; nearest cluster member {state.get('nearest', 'n/a')} SNPs; {state.get('members_within_ambiguous', 0)} members within {POLICY['ambiguous_max_snps']} SNPs."
def plan(state):                                           # the next tool, chosen from the current state (a language model makes this choice in ex1_7)
    if "qc" not in state: return "qc_check"
    if not state["qc"]["passes_quality_gate"]: return "stop_and_report" if "decision" not in state else None
    if "nearest" not in state: return "nearest_distance"
    if "decision" not in state: return "decide"
    if "note" not in state: return "draft_note"
    return None
def run_agent(name):
    state = dict(name=name); trace = []
    for step in range(1, POLICY["max_steps"] + 1):
        tool = plan(state)
        if tool is None: break
        if tool == "qc_check": obs = qc_check(name); state["qc"] = obs; state["callable"] = obs["callable_fraction"]
        elif tool == "stop_and_report": state["decision"] = DECISIONS[3]; state["needs_human"] = True; obs = {"decision": state["decision"]}
        elif tool == "nearest_distance": obs = nearest_distance(name); state.update(nearest=obs["nearest_snps"], members_within_link=obs["cluster_members_within_2_snps"], members_within_ambiguous=obs["cluster_members_within_5_snps"])
        elif tool == "decide": state["decision"] = decision_by_policy(state["nearest"]); state["needs_human"] = state["decision"] == DECISIONS[1]; obs = {"decision": state["decision"]}
        elif tool == "draft_note": state["note"] = draft_note(state, state["decision"]); obs = {"note": state["note"]}
        trace.append(dict(step=step, tool=tool, observation=json.dumps(obs, default=str)))
        if tool == "stop_and_report": break
    return state, trace
if __name__ == "__main__":
    rows = []; lines = []
    for name, g, truth, c in CASES:
        st, tr = run_agent(name); lines.append(f"== {name} (truth: {truth}, callable {c:.2f})"); lines += [f"   step {t['step']}: {t['tool']:<17} -> {t['observation']}" for t in tr]
        rows.append(dict(case=name, truth=truth, callable=c, nearest_snps=st.get("nearest", np.nan), steps=len(tr), tools=" > ".join(t["tool"] for t in tr), decision=st["decision"], needs_human=st["needs_human"]))
    R = pd.DataFrame(rows); R.to_csv("table_agent_trace.csv", index=False); open("agent_trace.txt", "w").write("\n".join(lines)); print(R[["case", "truth", "callable", "nearest_snps", "steps", "decision", "needs_human"]].to_string(index=False))
