#!/usr/bin/env python3
"""Report which packages and tools are installed and which chapters are ready (used by scripts/check-install.sh).
   Run it with the Python of the environment you want to check:  .venv/bin/python tools/check_install.py --chapters 12,13"""
import argparse, importlib.metadata as md, os, platform, shutil, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import needed as N
ap = argparse.ArgumentParser(); ap.add_argument("--chapters", default="1-13"); A = ap.parse_args()
spec = N.load(); chs = N.chapter_list(A.chapters); need = N.needed(chs, spec); pins = N.pins(); rows = []; missing = set()
v = platform.python_version(); rows.append(("python", "OK" if tuple(map(int, v.split(".")[:2])) >= (3, 12) else "OLD", v, "3.12.3"))
for p in need["python"]:
    try: f = md.version(p)
    except md.PackageNotFoundError: f = ""
    t = pins.get(p.lower(), ""); st = "MISSING" if not f else ("DIFFERENT" if t and f != t else "OK")
    if st == "MISSING": missing.add(p)
    rows.append((f"py: {p}", st, f, t))
rs = shutil.which("Rscript")
if need["r"]:
    if not rs:
        rows.append(("R", "MISSING", "", "4.3.3")); missing.update(need["r"])
    else:
        def r_eval(expr):
            try: return subprocess.run([rs, "-e", expr], capture_output=True, text=True, timeout=120).stdout.strip()
            except Exception: return ""
        rows.append(("R", "OK", r_eval("cat(as.character(getRversion()))"), "4.3.3"))
        for p in need["r"]:
            f = r_eval(f"cat(tryCatch(as.character(packageVersion('{p}')), error=function(e) ''))")
            if not f: missing.add(p)
            rows.append((f"R: {p}", "OK" if f else "MISSING", f, ""))
for t in dict.fromkeys(["pandoc"] + need["quick"] + need["rebuild"]):
    names = ["plink1.9", "plink"] if t == "plink1.9" else [t]; found = next((shutil.which(n) for n in names if shutil.which(n)), None)
    kind = "needed" if t in need["quick"] else ("optional, for --docx" if t == "pandoc" else "rebuild only")
    if not found and kind == "needed": missing.add(t)
    rows.append((f"tool: {t} ({kind})", "OK" if found else ("MISSING" if kind == "needed" else "not installed"), found or "", ""))
w0 = max(len(r[0]) for r in rows) + 2; w1 = max(len(r[1]) for r in rows) + 2; w2 = max(len(r[2]) for r in rows) + 2
print(f"{'Item':<{w0}}{'Status':<{w1}}{'Found':<{w2}}Tested"); print("-" * (w0 + w1 + w2 + 8))
for r in rows: print(f"{r[0]:<{w0}}{r[1]:<{w1}}{r[2]:<{w2}}{r[3]}")
print(f"\nTested on: {spec['tested_on']}"); bad = 0
for c in chs:
    e = spec["chapters"][str(c)]; reqs = list(e["python"]) + list(spec["base"]["python"]) + list(e["r"]) + list(e["quick_tools"]) + []; lack = sorted({x for x in reqs if x in missing})
    if lack: print(f"Chapter {c:>2}: NOT READY, missing {', '.join(lack)}"); bad += 1
    else: print(f"Chapter {c:>2}: ready for the quick check ({e['title']})")
sys.exit(1 if bad else 0)
