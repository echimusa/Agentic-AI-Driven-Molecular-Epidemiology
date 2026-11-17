"""
Chapter 12, Hands-on 9 - Check a surveillance protocol: the arithmetic and the completeness that a script can verify.
Two protocols for the same aim, written as data. Protocol A is consistent. Protocol B carries six planted faults (the negative control of the check); one of them, a site that is allocated too much, breaks two rules. Ten rules test the numbers and the presence of the decisions;
none of them judges whether the aim, the consent route or the sharing terms are right, which stay with people. The rules use the results of Hands-on 3, 5 and 8.
Input: table_privacy.csv, table_embargo.csv. Output: table_protocol_check.csv.
"""
import math
import pandas as pd
detect_n = lambda p, c: math.ceil(math.log(1 - c) / math.log(1 - p))
AIMS = {"detect new variants", "estimate prevalence", "investigate outbreaks", "monitor antimicrobial resistance"}; STANDARDS = {"PHA4GE contextual data specification", "INSDC pathogen package"}
def protocol(**kw):
    p = dict(name="", aim="detect new variants", target_frequency=0.02, confidence=0.95, weekly_capacity=96, sites={"North": 0.40, "Centre": 0.35, "South": 0.25}, population_share={"North": 0.38, "Centre": 0.37, "South": 0.25},
             weeks_to_detect_stated=2, turnaround_target_days=14, routine_turnaround_target_days=14, metadata=["collection date", "place", "specimen type", "host"], metadata_standard="PHA4GE contextual data specification",
             release_tier="open", release_date_resolution="month", release_place_resolution="country", consent_basis="public-health surveillance mandate with ethics review", credit_rule="submitting laboratory named as author on outputs", embargo_days=30, max_embargo_days=30, review_date="2027-06-30")
    p.update(kw); return p
A = protocol(name="A: consistent")
B = protocol(name="B: six planted faults", weekly_capacity=96, sites={"North": 0.50, "Centre": 0.35, "South": 0.25}, weeks_to_detect_stated=1, metadata=["collection date", "place"], release_date_resolution="day", release_place_resolution="district", embargo_days=120, review_date=None)
priv = pd.read_csv("table_privacy.csv"); emb = pd.read_csv("table_embargo.csv").set_index("embargo_days")
day_unique = float(priv[(priv.records == "day-level dates") & (priv.qi_set == "place + date + serovar") & (priv.resolution == "day")].unique.iloc[0])
def check(p):
    r = []; add = lambda i, rule, ok, detail: r.append(dict(protocol=p["name"], rule=i, check=rule, result="pass" if ok else "FAIL", detail=detail))
    n = detect_n(p["target_frequency"], p["confidence"]); weeks = n / p["weekly_capacity"]
    add("R1", "The aim is one of the four aims of the Session 12 table", p["aim"] in AIMS, p["aim"])
    add("R2", "The target frequency and the confidence are stated, and the sequences needed follow from them", 0 < p["target_frequency"] < 1 and 0 < p["confidence"] < 1, f"{n} sequences for {p['target_frequency']:.0%} at {p['confidence']:.0%}")
    add("R3", "The weeks to detect stated in the protocol equal the sequences needed divided by the weekly capacity, rounded up", p["weeks_to_detect_stated"] == math.ceil(weeks), f"computed {weeks:.2f} weeks, so {math.ceil(weeks)}; stated {p['weeks_to_detect_stated']}")
    tot = sum(p["sites"].values()); add("R4", "The capacity allocated to the sites adds up to 100%", abs(tot - 1) < 1e-9, f"sum {tot:.2f}")
    dev = max(abs(p["sites"][s] - p["population_share"][s]) for s in p["sites"]); add("R5", "For a representative aim no site departs from its population share by more than 10 percentage points", dev <= 0.10, f"largest departure {dev * 100:.0f} points")
    add("R6", "Turnaround targets are stated, and the outbreak target is not slower than the routine one", p["turnaround_target_days"] is not None and p["turnaround_target_days"] <= p["routine_turnaround_target_days"], f"{p['turnaround_target_days']} against {p['routine_turnaround_target_days']} days")
    need = {"collection date", "place", "specimen type", "host"}; miss = need - set(p["metadata"]); add("R7", "The minimum metadata names date, place, specimen type and host, and a standard is named", not miss and p["metadata_standard"] in STANDARDS, f"missing: {sorted(miss) or 'none'}; standard: {p['metadata_standard']}")
    risky = p["release_tier"] == "open" and p["release_date_resolution"] == "day" and p["release_place_resolution"] != "country"; add("R8", "An open release does not combine day-level dates with place below country level (Hands-on 3)", not risky, f"{p['release_tier']} release, {p['release_date_resolution']} dates, place at {p['release_place_resolution']} level; {day_unique:.1f}% of day-level human records were unique on place, date and serovar")
    add("R9", "The embargo is within the maximum the protocol allows, and its delay is stated", p["embargo_days"] <= p["max_embargo_days"], f"{p['embargo_days']} days against {p['max_embargo_days']}; Hands-on 8 median delay with the closest tested embargo ({min(emb.index, key=lambda e: abs(e - p['embargo_days']))} days): {emb.median_delay_days.loc[min(emb.index, key=lambda e: abs(e - p['embargo_days']))]:.0f} days")
    add("R10", "A consent or legal basis, a credit rule and a review date are present", bool(p["consent_basis"]) and bool(p["credit_rule"]) and bool(p["review_date"]), f"consent: {bool(p['consent_basis'])}; credit: {bool(p['credit_rule'])}; review: {p['review_date']}")
    return r
T = pd.DataFrame(check(A) + check(B)); T.to_csv("table_protocol_check.csv", index=False)
print(T[["protocol", "rule", "result", "detail"]].to_string(index=False, max_colwidth=100)); print(T.groupby("protocol").result.apply(lambda s: f"{(s == 'pass').sum()} of {len(s)} rules pass").to_string())
