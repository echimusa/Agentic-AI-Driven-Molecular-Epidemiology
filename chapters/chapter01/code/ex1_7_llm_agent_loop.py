"""
Chapter 1, agentic AI foundations - the cluster-triage agent with a language model choosing the tools.
Same scenario, tools and policy as ex1_6_agent_loop.py (triage_core.py). What changes is the planner: instead of a rule table, a language model sees the task, calls tools through the Anthropic Messages API tool-use
interface, reads the results and decides when to finish. What does not change is everything that carries risk, which stays in code and is enforced whatever the model says:
  G1  a tool call must name a known tool and the case under triage;
  G2  nearest_distance is blocked until qc_check has passed the quality gate;
  G3  the final decision must agree with the policy for the numbers the tools returned; two rejected decisions send the case to a person;
  G4  at most 8 tool calls per case; more send the case to a person;
  G5  a model that stops without calling finish sends the case to a person;
  G6  one tool call at a time (parallel calls disabled), and every call, result and guardrail event is logged.
USAGE
  export ANTHROPIC_API_KEY=...            # your key; it is read from the environment and never printed
  pip install anthropic pandas numpy
  python ex1_7_llm_agent_loop.py                         # live run on all nine cases (about 40 model calls)
  python ex1_7_llm_agent_loop.py --model claude-haiku-5-5 --runs 3
  python ex1_7_llm_agent_loop.py --mock                  # no key: a scripted stand-in replays well-behaved and badly-behaved 'models' to test the guardrails
The model name defaults to claude-sonnet-5-5 (override with --model or CLAUDE_MODEL). Model outputs vary between runs and models; use --runs to repeat and compare.
Output: table_llm_agent_trace.csv, llm_agent_trace.jsonl (live runs) or table_mock_agent_trace.csv, mock_agent_trace.jsonl (--mock).
"""
import argparse, json, os, sys
from types import SimpleNamespace as NS
import numpy as np, pandas as pd
from triage_core import *
SYSTEM = f"""You triage newly sequenced viral genomes for a public health team. For the case you are given, decide whether the new genome belongs to a known outbreak cluster of {len(CLUSTER)} sequenced cases.
Work with the tools only, and use only numbers that tools return; never guess a number. Call qc_check first. Call nearest_distance only if the quality gate passed.
The team's policy, which you cannot change: a nearest cluster member at most {POLICY['link_max_snps']} SNPs away means 'consistent with the cluster'; {POLICY['link_max_snps'] + 1} to {POLICY['ambiguous_max_snps']} SNPs means 'ambiguous: escalate to a person'; more than {POLICY['ambiguous_max_snps']} SNPs means 'not linked to the cluster'; a callable fraction below {POLICY['min_callable']} means the genome is not assessed and must be resequenced.
Finish by calling finish exactly once, with one of the allowed decisions and a one-sentence rationale that cites the numbers the tools returned. If anything is unclear or a tool is blocked, choose the escalation decision."""
TOOLS = [
    {"name": "qc_check", "description": "Quality gate for the new genome: returns the callable fraction (share of genome positions with a usable base call) and whether it passes the gate.", "input_schema": {"type": "object", "properties": {"case_id": {"type": "string"}}, "required": ["case_id"]}},
    {"name": "nearest_distance", "description": "SNP distance from the new genome to its nearest member of the known cluster, and how many cluster members are within 2 and within 5 SNPs. Only allowed after the quality gate has passed.", "input_schema": {"type": "object", "properties": {"case_id": {"type": "string"}}, "required": ["case_id"]}},
    {"name": "finish", "description": "Record the final decision for the case and end the task. Call exactly once.", "input_schema": {"type": "object", "properties": {"case_id": {"type": "string"}, "decision": {"type": "string", "enum": DECISIONS}, "rationale": {"type": "string"}}, "required": ["case_id", "decision", "rationale"]}}]
# ------------------------------------------------------------------ the guarded loop
def run_case(client, model, name, max_tokens=600):
    st = dict(qc=None, nearest=None, rejections=0, final=None, proposed=None, escalated_by=None); events = []; calls = 0; usage = [0, 0]
    messages = [{"role": "user", "content": f"Case {name}: triage this new genome against the known cluster."}]
    def handle(tool, inp):
        nonlocal calls
        if tool not in ("qc_check", "nearest_distance", "finish"): return {"error": f"unknown tool {tool}"}, True, "G1"
        if inp.get("case_id") != name: return {"error": f"case_id must be {name}"}, True, "G1"
        if tool == "qc_check": st["qc"] = qc_check(name); return st["qc"], False, None
        if tool == "nearest_distance":
            if st["qc"] is None or not st["qc"]["passes_quality_gate"]: return {"error": "blocked by policy: run qc_check and pass the quality gate first"}, True, "G2"
            st["nearest"] = nearest_distance(name); return st["nearest"], False, None
        dec = inp.get("decision"); st["proposed"] = dec                                                             # tool == "finish"
        if st["qc"] is not None and not st["qc"]["passes_quality_gate"]: expected = DECISIONS[3]
        elif st["nearest"] is not None: expected = decision_by_policy(st["nearest"]["nearest_snps"])
        else: expected = None
        if expected is None or dec != expected:
            st["rejections"] += 1; msg = "rejected by policy: the decision does not follow from the tool results obtained so far" + (" (run the missing tools)" if expected is None else "")
            if st["rejections"] >= 2: st["final"] = DECISIONS[1]; st["escalated_by"] = "G3"; return {"error": msg + "; escalated to a person"}, True, "G3"
            return {"error": msg}, True, "G3"
        st["final"] = dec; return {"recorded": dec}, False, None
    for turn in range(POLICY["max_steps"] + 2):
        r = client.messages.create(model=model, max_tokens=max_tokens, system=SYSTEM, tools=TOOLS, tool_choice={"type": "auto", "disable_parallel_tool_use": True}, messages=messages)
        if getattr(r, "usage", None): usage[0] += getattr(r.usage, "input_tokens", 0); usage[1] += getattr(r.usage, "output_tokens", 0)
        messages.append({"role": "assistant", "content": r.content})
        if r.stop_reason != "tool_use": st["final"] = st["final"] or DECISIONS[1]; st["escalated_by"] = st["escalated_by"] or "G5"; events.append(dict(event="stopped without finish", guardrail="G5")); break
        results = []
        for b in [b for b in r.content if b.type == "tool_use"]:
            calls += 1
            if calls > POLICY["max_steps"]: st["final"] = DECISIONS[1]; st["escalated_by"] = "G4"; calls = POLICY["max_steps"]; events.append(dict(event="step budget exceeded", guardrail="G4")); break
            res, err, g = handle(b.name, dict(b.input)); events.append(dict(step=calls, tool=b.name, input=dict(b.input), result=res, guardrail=g)); results.append({"type": "tool_result", "tool_use_id": b.id, "content": json.dumps(res), "is_error": err})
        if st["escalated_by"] or st["final"]: break
        messages.append({"role": "user", "content": results})
    if st["final"] is None: st["final"] = DECISIONS[1]; st["escalated_by"] = "G4"
    needs_human = st["final"] in (DECISIONS[1], DECISIONS[3]) or st["escalated_by"] is not None
    return dict(case=name, final=st["final"], needs_human=needs_human, tool_calls=calls, guardrail_events=";".join(sorted({e["guardrail"] for e in events if e.get("guardrail")})), escalated_by=st["escalated_by"] or "", model_proposed=st["proposed"] or "", input_tokens=usage[0], output_tokens=usage[1]), events
# ------------------------------------------------------------------ scripted stand-in for the model (--mock): tests the guardrails, says nothing about any real model
class ScriptedClient:
    def __init__(self):
        P = lambda c: decision_by_policy(nearest_distance(c)["nearest_snps"]) if qc_check(c)["passes_quality_gate"] else DECISIONS[3]
        T = lambda n, i: ("tool", n, i); F = lambda c, d: ("tool", "finish", dict(case_id=c, decision=d, rationale="scripted"))
        self.scripts = {}
        for c, *_ in CASES:
            q = T("qc_check", dict(case_id=c)); n = T("nearest_distance", dict(case_id=c))
            self.scripts[c] = [q, F(c, DECISIONS[3])] if not qc_check(c)["passes_quality_gate"] else [q, n, F(c, P(c))]                      # a well-behaved model
        self.scripts["new-2"] = [T("nearest_distance", dict(case_id="new-2")), T("qc_check", dict(case_id="new-2")), T("nearest_distance", dict(case_id="new-2")), F("new-2", P("new-2"))]       # tries the distance before the gate (G2)
        self.scripts["new-3"] = [T("qc_check", dict(case_id="new-3")), T("nearest_distance", dict(case_id="new-3")), F("new-3", DECISIONS[2]), F("new-3", P("new-3"))]                           # one wrong decision, then right (G3)
        self.scripts["new-4"] = [T("qc_check", dict(case_id="new-4")), T("nearest_distance", dict(case_id="new-4")), F("new-4", DECISIONS[0]), F("new-4", DECISIONS[0])]                          # two wrong decisions: escalated (G3)
        self.scripts["new-5"] = [T("qc_check", dict(case_id="new-5")), T("nearest_distance", dict(case_id="new-5")), ("text", "The genome looks closely related to the cluster.")]                       # stops without finish (G5)
        self.scripts["new-6"] = [T("qc_check", dict(case_id="new-6"))] * 9                                                                                                                         # loops on one tool: step budget (G4)
        self.scripts["new-8"] = [T("qc_check", dict(case_id="new-9")), T("qc_check", dict(case_id="new-8")), T("nearest_distance", dict(case_id="new-8")), F("new-8", P("new-8"))]                  # asks about another case (G1)
        self.messages = self
    def create(self, **kw):
        msgs = kw["messages"]; name = msgs[0]["content"].split()[1].rstrip(":"); k = sum(1 for m in msgs if m["role"] == "assistant"); steps = self.scripts[name]; step = steps[min(k, len(steps) - 1)] if k < len(steps) else ("text", "done")
        usage = NS(input_tokens=0, output_tokens=0)
        if step[0] == "text": return NS(stop_reason="end_turn", content=[NS(type="text", text=step[1])], usage=usage)
        return NS(stop_reason="tool_use", content=[NS(type="tool_use", id=f"toolu_mock_{name}_{k}", name=step[1], input=step[2])], usage=usage)
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--mock", action="store_true"); ap.add_argument("--model", default=os.environ.get("CLAUDE_MODEL", "claude-sonnet-5-5")); ap.add_argument("--runs", type=int, default=1); ap.add_argument("--cases", default="")
    a = ap.parse_args(); names = [c[0] for c in CASES if not a.cases or c[0] in a.cases.split(",")]
    if a.mock: client = ScriptedClient(); tag = "mock"
    else:
        if not os.environ.get("ANTHROPIC_API_KEY"): sys.exit("Set ANTHROPIC_API_KEY, or run with --mock to test the guardrails without a model.")
        try: import anthropic
        except ImportError: sys.exit("pip install anthropic")
        client = anthropic.Anthropic(); tag = "llm"
    rows = []; log = open(f"{'mock' if a.mock else 'llm'}_agent_trace.jsonl", "w")
    for run in range(a.runs):
        for name in names:
            row, events = run_case(client, a.model, name); truth = BY_NAME[name][1][2]; rule, rule_h, _ = rule_agent(name); row.update(run=run + 1, truth=truth, rule_agent_decision=rule, agrees_with_rule_agent=row["final"] == rule); rows.append(row)
            log.write(json.dumps(dict(run=run + 1, case=name, events=events, summary=row), default=str) + "\n")
    R = pd.DataFrame(rows); R.to_csv(f"table_{tag}_agent_trace.csv", index=False)
    print(R[["run", "case", "truth", "final", "needs_human", "tool_calls", "guardrail_events", "escalated_by", "agrees_with_rule_agent"]].to_string(index=False)); print(f"\nagreement with the rule-table agent: {R.agrees_with_rule_agent.mean():.0%}; sent to a person: {R.needs_human.mean():.0%}; guardrail events in {int((R.guardrail_events != '').sum())} of {len(R)} case runs")
if __name__ == "__main__": main()
