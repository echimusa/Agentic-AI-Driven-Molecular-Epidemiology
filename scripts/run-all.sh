#!/usr/bin/env bash
# Run the quick check (or the full rerun) of several chapters and print one summary table. Same job as run-all.ps1.
# Usage: scripts/run-all.sh [--chapters 1,2,4,12,13] [--full]
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
CH_ARG="1-13"; FULL=""
while [ $# -gt 0 ]; do case "$1" in --chapters) CH_ARG="${2:?--chapters needs a list}"; shift 2;; --full) FULL="--full"; shift;; -h|--help) sed -n 2,3p "$0"; exit 0;; *) die "unknown option $1";; esac; done
FAILS=0; SUMMARY=""
for c in $(parse_chapters "$CH_ARG"); do
  start=$(date +%s); "$ROOT/scripts/run-chapter.sh" "$c" $FULL; code=$?; secs=$(( $(date +%s) - start ))
  if [ $code -eq 0 ]; then res=PASS; else res=FAIL; FAILS=$((FAILS + 1)); fi
  SUMMARY+=$(printf '%-11s %-7s %s s\n' "Chapter $c" "$res" "$secs")$'\n'
done
printf '\n%-11s %-7s %s\n%s' "Chapter" "Result" "Time" "$SUMMARY"
[ $FAILS -eq 0 ]
