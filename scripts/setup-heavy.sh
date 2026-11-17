#!/usr/bin/env bash
# Download BEAST 2.7.7 (with the MODEL_SELECTION package) and AMRFinderPlus 4.2.7 for the heavy analyses of Chapters 6 and 9. Same job as setup-heavy.ps1.
# Usage: scripts/setup-heavy.sh [--skip-beast] [--skip-amrfinder] [--dry-run]    (needs Java, unzip and curl: scripts/setup.sh --heavy)
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
BEAST=1; AMR=1; DRY=0
while [ $# -gt 0 ]; do case "$1" in --skip-beast) BEAST=0; shift;; --skip-amrfinder) AMR=0; shift;; --dry-run) DRY=1; shift;; -h|--help) sed -n 2,3p "$0"; exit 0;; *) die "unknown option $1";; esac; done
run() { say "$*"; if [ "$DRY" = 1 ]; then echo "    (dry run: not executed)"; else eval "$*"; fi; }
BIN="$ROOT/tools/linux"; mkdir -p "$BIN/bin"
if [ "$BEAST" = 1 ]; then
  run "curl -L -o '$BIN/beast.tgz' https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz && tar -xzf '$BIN/beast.tgz' -C '$BIN' && rm '$BIN/beast.tgz'"
  run "for t in beast beauti treeannotator packagemanager loganalyser; do [ -f '$BIN/beast/bin/'\$t ] && ln -sf '$BIN/beast/bin/'\$t '$BIN/bin/'\$t; done; true"
  run "mkdir -p \"\$HOME/.beast/2.7/MODEL_SELECTION\" && curl -L -o '$BIN/ms.zip' https://github.com/BEAST2-Dev/model-selection/releases/download/v1.6.0/MODEL_SELECTION.addon.v1.6.2.zip && unzip -o '$BIN/ms.zip' -d \"\$HOME/.beast/2.7/MODEL_SELECTION\" >/dev/null && rm '$BIN/ms.zip'"
fi
if [ "$AMR" = 1 ]; then
  run "mkdir -p '$ROOT/amrf' && curl -L -o '$ROOT/amrf/amrf.tgz' https://github.com/ncbi/amr/releases/download/amrfinder_v4.2.7/amrfinder_binaries_v4.2.7.tar.gz && tar -xzf '$ROOT/amrf/amrf.tgz' -C '$ROOT/amrf' && rm '$ROOT/amrf/amrf.tgz'"
  echo "Then fetch the database files (docs/HEAVY-ANALYSES.md lists release 2026-10-29.1) into amrf/db and run: amrf/amrfinder_index amrf/db"
fi
echo "BEAST 2 programs are in $BIN/bin; scripts/_common.sh adds that folder to PATH for the other scripts."
