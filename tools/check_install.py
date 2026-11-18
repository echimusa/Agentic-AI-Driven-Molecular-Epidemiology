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
# heavy tools: only needed to repeat the BEAST 2 and AMRFinderPlus analyses of Chapters 6 and 9, so they never make a chapter "not ready"
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); want = {t for c in chs for t in spec["chapters"][str(c)].get("heavy_tools", [])}; heavy_todo = []
def first_line(cmd):
    try: out = subprocess.run(cmd, capture_output=True, text=True, timeout=60); t = (out.stdout + out.stderr).strip(); return t.splitlines()[0] if t else ""
    except Exception: return ""
def heavy_row(name, found, detail, tested, ok):
    rows.append((f"heavy: {name}", "OK" if ok else "not installed", detail if ok else "", tested))
    if not ok: heavy_todo.append(name)
if "beast2" in want:
    b = next((q for q in (os.path.join(root, "tools", "linux", "bin", "beast"), shutil.which("beast")) if q and os.path.exists(q)), None)
    heavy_row("BEAST 2 (Chapter 6)", b, b or "", "2.7.7", bool(b) and "2.7.7" in first_line([b, "-version"]))
    jre = os.path.join(root, "tools", "linux", "beast", "jre", "bin", "java"); jv = shutil.which("java")
    heavy_row("Java for BEAST 2", jre if os.path.exists(jre) else jv, "bundled with BEAST 2" if os.path.exists(jre) else (jv or ""), "17", os.path.exists(jre) or bool(jv))
    ms = os.path.join(os.path.expanduser("~"), ".beast", "2.7", "MODEL_SELECTION", "version.xml")
    heavy_row("BEAST 2 package MODEL_SELECTION", ms, os.path.dirname(ms), "1.6.2", os.path.exists(ms) and "1.6.2" in open(ms, errors="ignore").read(2000))
if "amrfinder" in want:
    a = os.path.join(root, "amrf", "amrfinder"); ver = first_line([a, "--version"]) if os.path.exists(a) else ""
    heavy_row("AMRFinderPlus (Chapter 9)", a, f"{a} {ver}", "4.2.7", ver == "4.2.7")
    vt = os.path.join(root, "amrf", "db", "version.txt"); dbv = open(vt).readline().strip() if os.path.exists(vt) else ""
    heavy_row("AMRFinderPlus database", vt, f"release {dbv}", "2026-10-29.1", dbv == "2026-10-29.1" and os.path.exists(os.path.join(root, "amrf", "db", ".indexed-2026-10-29.1")))
if "rgi" in want:
    r_ = shutil.which("rgi"); heavy_row("RGI (Chapter 9, optional, installed by hand)", r_, r_ or "", "6.0.8", bool(r_))
w0 = max(len(r[0]) for r in rows) + 2; w1 = max(len(r[1]) for r in rows) + 2; w2 = max(len(r[2]) for r in rows) + 2
print(f"{'Item':<{w0}}{'Status':<{w1}}{'Found':<{w2}}Tested"); print("-" * (w0 + w1 + w2 + 8))
for r in rows: print(f"{r[0]:<{w0}}{r[1]:<{w1}}{r[2]:<{w2}}{r[3]}")
print(f"\nTested on: {spec['tested_on']}"); bad = 0
for c in chs:
    e = spec["chapters"][str(c)]; reqs = list(e["python"]) + list(spec["base"]["python"]) + list(e["r"]) + list(e["quick_tools"]) + []; lack = sorted({x for x in reqs if x in missing})
    if lack: print(f"Chapter {c:>2}: NOT READY, missing {', '.join(lack)}"); bad += 1
    else: print(f"Chapter {c:>2}: ready for the quick check ({e['title']})")
todo = [h for h in heavy_todo if "RGI" not in h]
if todo: print("\nNot installed (needed only to repeat the BEAST 2 and AMRFinderPlus analyses): " + ", ".join(todo) + ".  Run scripts/setup-heavy.sh (or scripts/setup.sh --heavy) to check them and install what is missing.")
sys.exit(1 if bad else 0)
