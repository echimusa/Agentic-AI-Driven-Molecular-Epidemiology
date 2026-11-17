#!/usr/bin/env bash
# Report what is installed and which chapters are ready to run from the shipped data. Same job as check-install.ps1.
# Usage: scripts/check-install.sh [--chapters 12,13]
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
CH_ARG="1-13"
while [ $# -gt 0 ]; do case "$1" in --chapters) CH_ARG="${2:?--chapters needs a list}"; shift 2;; -h|--help) sed -n 2,3p "$0"; exit 0;; *) die "unknown option $1";; esac; done
set_env; PY="$(get_python)" || die "python3 not found; run scripts/setup.sh first"
exec "$PY" "$ROOT/tools/check_install.py" --chapters "$(parse_chapters "$CH_ARG" | tr ' ' ',')"
