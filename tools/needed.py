#!/usr/bin/env python3
"""Print what the chosen chapters need, as shell variable assignments (used by scripts/setup.sh and scripts/check-install.sh).
   python3 tools/needed.py 12,13 [--rebuild] [--heavy]"""
import argparse, json, os, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def chapter_list(text):
    out = set()
    for part in re.split(r"[,\s]+", text.strip()):
        if not part: continue
        m = re.fullmatch(r"(\d+)(?:-|\.\.)(\d+)", part)
        out.update(range(int(m.group(1)), int(m.group(2)) + 1) if m else [int(part)])
    bad = [c for c in out if not 1 <= c <= 13]
    if bad: sys.exit(f"chapters must be between 1 and 13, not {bad}")
    return sorted(out)
def load():
    return json.load(open(os.path.join(ROOT, "tools", "chapter_requirements.json")))
def pins():
    d = {}
    for line in open(os.path.join(ROOT, "requirements.txt"), encoding="utf8"):
        if "==" in line and not line.lstrip().startswith("#"): n, v = line.strip().split("=="); d[n.lower()] = v
    return d
def needed(chs, spec=None):
    spec = spec or load(); py = list(spec["base"]["python"]); r = []; quick = []; rebuild = []; heavy = []
    for c in chs:
        e = spec["chapters"][str(c)]
        py += [p for p in e["python"] if p not in py]; r += [p for p in e["r"] if p not in r]
        quick += [p for p in e["quick_tools"] if p not in quick]; rebuild += [p for p in e["rebuild_tools"] if p not in rebuild]; heavy += [p for p in e.get("heavy_tools", []) if p not in heavy]
    return dict(python=py, r=r, quick=quick, rebuild=rebuild, heavy=heavy)
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("chapters", nargs="?", default="1-13"); ap.add_argument("--rebuild", action="store_true"); ap.add_argument("--heavy", action="store_true"); A = ap.parse_args()
    spec = load(); chs = chapter_list(A.chapters); n = needed(chs, spec); pn = pins()
    apt = []
    for t in n["quick"] + (n["rebuild"] if A.rebuild else []) + (n["heavy"] if A.heavy else []):
        for p in spec["apt"].get(t, "").split():
            if p not in apt: apt.append(p)
    q = lambda xs: "'" + " ".join(xs) + "'"
    print(f"CHAPTERS={q(map(str, chs))}"); print(f"PY_NAMES={q(n['python'])}"); print(f"PY_SPECS={q([p + '==' + pn[p.lower()] if p.lower() in pn else p for p in n['python']])}")
    print(f"R_PKGS={q(n['r'])}"); print(f"QUICK_TOOLS={q(n['quick'])}"); print(f"REBUILD_TOOLS={q(n['rebuild'])}"); print(f"HEAVY_TOOLS={q(n['heavy'])}"); print(f"APT_EXTRA={q(apt)}")
