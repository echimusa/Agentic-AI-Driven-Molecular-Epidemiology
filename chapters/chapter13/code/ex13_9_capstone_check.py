"""
Chapter 13, Hands-on 9 - Check a capstone analysis plan and briefing: what a script can verify about the Session 14 workflow, deliverables and briefing.
Two plans for the same capstone outbreak, written as data. Plan A is consistent. Plan B carries six planted faults (the negative control of the check). Ten rules test form: the stages, the scripts, the written assumptions,
the uncertainty, the trees, the deliverable order, and the shape of the briefing. None of them judges whether the evidence supports the recommendation, which stays with the team and its tutors.
Output: table_capstone_check.csv.
"""
import pandas as pd
STAGES = ["QC and typing", "Clusters", "Trees and dates", "Transmission and spread", "Recommendations"]
def plan(**kw):
    p = dict(name="", stages={s: dict(script=f"scripts/{i + 1:02d}_{s.split()[0].lower()}.py", assumptions=[f"assumption for {s}"]) for i, s in enumerate(STAGES)}, trees="posterior sample of 100 trees",
             findings=[dict(text="Cluster of 9 isolates at 3 SNPs or fewer", interval="95% CrI 7-12", confidence="high"), dict(text="Ongoing transmission on ward 4", interval="probability 0.72", confidence="moderate")],
             deliverables=[("Analysis plan", 13), ("Team presentation", 15), ("Practical portfolio", 22), ("Report", 29)], first_sentence="Decision needed: whether to close ward 4 to admissions this week.", sections=["Decision", "Evidence", "Recommendation"], words=430, figures=1)
    p.update(kw); return p
A = plan(name="A: consistent")
B_stages = {s: dict(script=f"scripts/{i + 1:02d}_{s.split()[0].lower()}.py", assumptions=[f"assumption for {s}"]) for i, s in enumerate(STAGES)}; B_stages["Clusters"]["script"] = None
B = plan(name="B: six planted faults", stages=B_stages, trees="one maximum-likelihood tree", findings=[dict(text="Cluster of 9 isolates at 3 SNPs or fewer", interval=None, confidence="high"), dict(text="Ongoing transmission on ward 4", interval="probability 0.72", confidence="moderate")],
         deliverables=[("Analysis plan", 13), ("Team presentation", 15), ("Practical portfolio", 29), ("Report", 22)], first_sentence="This briefing summarises a genomic investigation of ward 4.", words=620, figures=2)
def check(p):
    r = []; add = lambda i, rule, ok, detail: r.append(dict(plan=p["name"], rule=i, check=rule, result="pass" if ok else "FAIL", detail=detail))
    add("C1", "All five stages of the capstone workflow are present", list(p["stages"]) == STAGES, f"{len(p['stages'])} stages")
    noscript = [s for s, v in p["stages"].items() if not v["script"]]; add("C2", "Every stage names a script in scripts/ (no manual clicks without records)", not noscript, f"no script for: {noscript or 'none'}")
    noass = [s for s, v in p["stages"].items() if not v["assumptions"]]; add("C3", "Every stage lists at least one written assumption", not noass, f"none for: {noass or 'none'}")
    nointerval = [f["text"][:30] for f in p["findings"] if not f["interval"]]; add("C4", "Every quantitative finding reports an interval or a probability", not nointerval, f"missing for: {nointerval or 'none'}")
    add("C5", "The trees are a sample or more than one, not a single best tree", "posterior" in p["trees"] or "sample" in p["trees"] or "bootstrap" in p["trees"], p["trees"])
    due = [d for _, d in p["deliverables"]]; add("C6", "Deliverables fall due in order: plan, presentation, portfolio, report", due == sorted(due), f"days {due}")
    add("C7", "The briefing leads with the decision needed", p["first_sentence"].lower().startswith("decision needed"), p["first_sentence"][:60])
    nocon = [f["text"][:30] for f in p["findings"] if f["confidence"] not in ("high", "moderate", "low")]; add("C8", "Every finding carries a confidence level", not nocon, f"missing for: {nocon or 'none'}")
    add("C9", "Evidence and recommendation are separate sections", "Evidence" in p["sections"] and "Recommendation" in p["sections"], ", ".join(p["sections"]))
    add("C10", "One page (at most 500 words) and one figure", p["words"] <= 500 and p["figures"] == 1, f"{p['words']} words, {p['figures']} figures")
    return r
T = pd.DataFrame(check(A) + check(B)); T.to_csv("table_capstone_check.csv", index=False)
print(T[["plan", "rule", "result", "detail"]].to_string(index=False, max_colwidth=70)); print(T.groupby("plan").result.apply(lambda s: f"{(s == 'pass').sum()} of {len(s)} rules pass").to_string())
