#!/usr/bin/env python3
"""Run the checks or one exercise of a chapter of the book, from the shipped data.

Examples (from the repository root, with the environment of scripts/setup.ps1 active):
    python tools/run_chapter.py 12                 # quick checks: recompute every number of the chapter (about a minute)
    python tools/run_chapter.py 12 --full          # also rerun every script and compare each result table byte for byte
    python tools/run_chapter.py 13 --script ex13_5_sir_gi.py   # run one exercise and keep its output in chapters/chapter13/.work
    python tools/run_chapter.py 6 --full --sim     # add the long TreeTime reruns of Chapter 6
    python tools/run_chapter.py 5 --doi            # also resolve every DOI in Crossref (needs the network)
Chapters 1 to 7 are run in a flat work folder (code and data side by side), as their notes describe; chapters 8 to 13 are run in review/.
"""
import argparse, glob, os, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("chapter", type=int, help="chapter number, 1 to 13")
ap.add_argument("--full", action="store_true", help="rerun every script and compare the result tables")
ap.add_argument("--doi", action="store_true", help="resolve every DOI of the chapter in Crossref (network)")
ap.add_argument("--sim", action="store_true", help="include the long simulations or reruns (chapters 5 to 8)")
ap.add_argument("--script", help="run one exercise script (file name in code/) in the work folder instead of the validator")
ap.add_argument("--offline", action="store_true", help="block the network for the run (proves it works from the shipped data)")
ap.add_argument("--docx", help="chapter Word file to check (default: the one in manuscript/)")
A = ap.parse_args()
if not 1 <= A.chapter <= 13: sys.exit("chapter must be between 1 and 13")
NN = f"{A.chapter:02d}"; CH = os.path.join(ROOT, "chapters", f"chapter{NN}"); n = A.chapter
if not os.path.isdir(CH): sys.exit(f"missing folder {CH}")
docx = A.docx or next(iter(sorted(glob.glob(os.path.join(CH, "manuscript", "Chapter_*.docx")))), None)
env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
if A.offline: env.update({k: "http://127.0.0.1:9" for k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY")})

def flat_workdir():
    work = os.path.join(CH, ".work"); shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    for sub in ("code", "data"):
        for f in glob.glob(os.path.join(CH, sub, "*")):
            if os.path.isfile(f): shutil.copy(f, work)
    for f in glob.glob(os.path.join(CH, "review", "validate_ch*.py")): shutil.copy(f, work)
    # the scripts write their figures into folders such as figs13/ and expect them to exist
    import re
    for f in glob.glob(os.path.join(CH, "code", "*.py")) + glob.glob(os.path.join(CH, "code", "*.R")):
        for m in re.finditer(r"""["']([A-Za-z_]*figs?[A-Za-z_0-9]*)/""", open(f, encoding="utf8", errors="ignore").read()): os.makedirs(os.path.join(work, m.group(1)), exist_ok=True)
    return work

if A.script:
    path = os.path.join(CH, "code", A.script)
    if not os.path.isfile(path): sys.exit(f"no such script: {path}\navailable: {', '.join(sorted(os.path.basename(f) for f in glob.glob(os.path.join(CH, 'code', 'ex*'))))}")
    work = flat_workdir()
    if os.path.isdir(os.path.join(ROOT, "amrf")):  # AMRFinderPlus installed by scripts/setup-heavy.sh: the exercises of Chapter 9 call ./amrf/amrfinder and amrf/db
        try: os.symlink(os.path.join(ROOT, "amrf"), os.path.join(work, "amrf"))
        except OSError: pass
    cmd = ["Rscript", A.script] if A.script.lower().endswith(".r") else [sys.executable, A.script]
    before = {f: os.path.getmtime(f) for f in glob.glob(os.path.join(work, "**", "*"), recursive=True) if os.path.isfile(f)}
    print(f"[chapter {n}] running {A.script} in {work}"); rc = subprocess.call(cmd, cwd=work, env=env)
    new = [os.path.relpath(f, work) for f in glob.glob(os.path.join(work, "**", "*"), recursive=True) if os.path.isfile(f) and (f not in before or os.path.getmtime(f) != before[f])]
    print(f"[chapter {n}] finished with exit code {rc}; {len(new)} file(s) written or changed in {work}: {', '.join(sorted(new)[:12])}{' ...' if len(new) > 12 else ''}"); sys.exit(rc)

val = os.path.join(CH, "review", f"validate_ch{n}.py")
if not os.path.isfile(val): sys.exit(f"missing validator {val}")
rec = os.path.join(CH, "review")
# A run without the chapter's Word file cannot produce the chapter-text rows of the committed record in review/, so it must not overwrite that record
# (it would leave the repository modified for everyone who runs a check): keep the record and put this run's table in .work instead.
keep = {} if docx else {f: open(f, "rb").read() for f in glob.glob(os.path.join(rec, "validation_results*.csv"))}
args = ["--chapter", docx] if docx else []
if not docx: print(f"[chapter {n}] no Word file given: the checks of the chapter text are skipped (the Word files are not in the repository; use --docx FILE)")
for flag, on in (("--full", A.full), ("--doi", A.doi), ("--sim", A.sim)):
    if on: args.append(flag)
if n <= 7:
    cwd = flat_workdir(); target = os.path.join(cwd, f"validate_ch{n}.py")
else:
    cwd = os.path.join(CH, "review"); target = val
print(f"[chapter {n}] {' '.join(['python', os.path.basename(target)] + args)}   (in {cwd})")
rc = subprocess.call([sys.executable, target] + args, cwd=cwd, env=env)
if n <= 7:
    for f in glob.glob(os.path.join(cwd, "validation_results*.csv")): shutil.copy(f, os.path.join(CH, "review", os.path.basename(f)))
    if rc == 0: shutil.rmtree(cwd, ignore_errors=True)
if keep:
    out = os.path.join(CH, ".work"); os.makedirs(out, exist_ok=True)
    for f, old in keep.items():
        shutil.copy(f, os.path.join(out, os.path.basename(f)))
        with open(f, "wb") as h: h.write(old)
    print(f"[chapter {n}] the table of this run is in {out}; the record in review/ (which includes the checks of the chapter text) was left unchanged. Use --docx FILE to refresh the record.")
print(f"[chapter {n}] exit code {rc} ({'all checks passed' if rc == 0 else 'see the lines marked FAIL above'})"); sys.exit(rc)
