#!/usr/bin/env bash
# Install what the exercises of the chosen chapters need on Ubuntu (or WSL2 Ubuntu): system packages, a Python environment in .venv, and R packages in R/library.
# Same job as setup.ps1.
usage() { cat <<'USAGE'
Usage: scripts/setup.sh [--chapters 12,13] [--rebuild] [--heavy] [--skip-system] [--here] [--dry-run]
  --chapters LIST   chapters to prepare, e.g. 12,13 or 1-5 (default 1-13)
  --rebuild         also install MAFFT, FastTree, bcftools, tabix and PLINK (to rebuild inputs from public sources)
  --heavy           also install Java, BLAST and HMMER for BEAST 2 and AMRFinderPlus (then run scripts/setup-heavy.sh)
  --skip-system     do not use apt (Python 3.12, R, pandoc and the tools are already installed)
  --here            WSL2 only: set up in this folder even if it is on a Windows drive (by default it is copied to the Linux home folder first)
  --dry-run         print every step and run none
USAGE
}
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
CHAPTERS_ARG="1-13"; REBUILD=""; HEAVY=""; SKIP=0; DRY=0; HERE=0; ORIG=("$@")
while [ $# -gt 0 ]; do case "$1" in
  --chapters) CHAPTERS_ARG="${2:?--chapters needs a list}"; shift 2;;
  --rebuild) REBUILD="--rebuild"; shift;; --heavy) HEAVY="--heavy"; shift;;
  --skip-system) SKIP=1; shift;; --here) HERE=1; shift;; --dry-run) DRY=1; shift;;
  -h|--help) usage; exit 0;; *) usage; die "unknown option $1";; esac; done
# WSL2 on a Windows drive (/mnt/c/...): slow, and files saved on Windows can have CRLF line endings. Copy the repository to the Linux home folder and run there.
if grep -qi microsoft /proc/version 2>/dev/null && [[ "$ROOT" == /mnt/* ]]; then
  if [ "$HERE" = 0 ]; then
    say "this copy is on a Windows drive ($ROOT): copying it to the Linux home folder first (use --here to stay)"
    MEPI_SRC="$ROOT" exec bash <(tr -d '\r' < "$ROOT/scripts/wsl-handoff.sh") "${ORIG[@]}"
  else warn "running on a Windows drive: this is slow, and files saved on Windows can have CRLF line endings"; fi
fi
command -v python3 >/dev/null 2>&1 || [ "$SKIP" = 0 ] || die "python3 not found"
CH="$(parse_chapters "$CHAPTERS_ARG")"
if command -v python3 >/dev/null 2>&1; then eval "$(python3 "$ROOT/tools/needed.py" "${CH// /,}" $REBUILD $HEAVY)"
else # python3 is missing: install it first, then read the requirements
  say "python3 is not installed; installing it first"; [ "$DRY" = 1 ] || { sudo apt-get update && sudo apt-get install -y python3 python3-venv python3-pip; }
  eval "$(python3 "$ROOT/tools/needed.py" "${CH// /,}" $REBUILD $HEAVY)"; fi
run() { say "$*"; if [ "$DRY" = 1 ]; then echo "    (dry run: not executed)"; else eval "$*"; fi; }
echo "Chapters: $CHAPTERS"; echo "Python packages: $PY_NAMES"; echo "R packages: ${R_PKGS:-(none)}"; echo "Tools needed for the quick check: ${QUICK_TOOLS:-(none)}"
[ "$(uname -s)" = "Linux" ] || die "these scripts are for Ubuntu (or WSL2 Ubuntu); on Windows see docs/INSTALL-Windows.md"
if [ "$SKIP" = 0 ]; then
  SUDO=""; [ "$(id -u)" = 0 ] || SUDO="sudo "
  APT="python3 python3-venv python3-pip python3-dev build-essential zlib1g-dev r-base r-base-dev pandoc curl unzip git libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libjpeg-dev $APT_EXTRA"
  run "${SUDO}apt-get update"
  run "${SUDO}env DEBIAN_FRONTEND=noninteractive apt-get install -y $APT"
fi
run "[ -x '$VENV_PY' ] || python3 -m venv '$ROOT/.venv'"
run "'$VENV_PY' -m pip install --upgrade pip"
FAILED=""
for spec in $PY_SPECS; do
  say "pip install $spec"
  if [ "$DRY" = 1 ]; then echo "    (dry run: not executed)"; elif ! "$VENV_PY" -m pip install "$spec"; then FAILED="$FAILED $spec"; fi
done
[ -z "$FAILED" ] || warn "pip could not install:$FAILED (see docs/TROUBLESHOOTING.md)"
if [ -n "$R_PKGS" ]; then
  say "install R packages: $R_PKGS"
  if [ "$DRY" = 1 ]; then echo "    (dry run: not executed)"; else set_env; command -v Rscript >/dev/null 2>&1 || die "Rscript not found; install r-base"; Rscript "$ROOT/r-packages.R" $R_PKGS || warn "some R packages did not install; see the lines marked MISSING and docs/TROUBLESHOOTING.md"; fi
fi
if [ -n "$HEAVY" ]; then  # BEAST 2, Java, AMRFinderPlus and its database: check each one, install what is missing
  say "check BEAST 2, Java and AMRFinderPlus, and install what is missing"
  if [ "$DRY" = 1 ]; then bash "$ROOT/scripts/setup-heavy.sh" --chapters "${CH// /,}" --dry-run || true
  else bash "$ROOT/scripts/setup-heavy.sh" --chapters "${CH// /,}" || warn "BEAST 2 or AMRFinderPlus is not fully installed: read the messages above, fix the cause and run scripts/setup-heavy.sh again"; fi
fi
if [ "$DRY" = 0 ]; then say "check the installation"; "$ROOT/scripts/check-install.sh" --chapters "${CH// /,}" || true; fi
