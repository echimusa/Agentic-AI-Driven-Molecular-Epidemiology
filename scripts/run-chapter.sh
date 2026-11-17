#!/usr/bin/env bash
# Check one chapter, or run one exercise, from the shipped data. Same job as run-chapter.ps1.
usage() { cat <<'USAGE'
Usage: scripts/run-chapter.sh CHAPTER [--full] [--doi] [--sim] [--offline] [--script FILE] [--docx FILE]
  scripts/run-chapter.sh 12                         recompute every number of the chapter and compare it with the text
  scripts/run-chapter.sh 12 --full                  also rerun every script and compare each result table byte for byte
  scripts/run-chapter.sh 13 --script ex13_5_sir_gi.py   run one exercise; output stays in chapters/chapter13/.work
  scripts/run-chapter.sh 6 --full --sim             include the long TreeTime reruns of Chapter 6
  --doi      resolve every DOI of the chapter in Crossref (network)      --offline  block the network for the run
  --docx FILE   also compare the text of your copy of the chapter's Word file with the results (the Word files are not in the repository)
USAGE
}
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
[ $# -ge 1 ] || { usage; exit 2; }
case "$1" in -h|--help) usage; exit 0;; esac
CH="$1"; shift; [[ $CH =~ ^[0-9]+$ ]] && (( CH >= 1 && CH <= 13 )) || die "the chapter must be a number from 1 to 13"
ARGS=()
while [ $# -gt 0 ]; do case "$1" in
  --full|--doi|--sim|--offline) ARGS+=("$1"); shift;;
  --script) ARGS+=(--script "${2:?--script needs a file name}"); shift 2;;
  --docx) ARGS+=(--docx "${2:?--docx needs a file name}"); shift 2;;
  -h|--help) usage; exit 0;; *) usage; die "unknown option $1";; esac; done
set_env; PY="$(get_python)" || die "python3 not found; run scripts/setup.sh first"
cd "$ROOT"; exec "$PY" tools/run_chapter.py "$CH" "${ARGS[@]}"
