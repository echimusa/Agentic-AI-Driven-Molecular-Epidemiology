#!/usr/bin/env bash
# Shared helpers for the Bash scripts (Ubuntu, and WSL2 Ubuntu). Source this file; do not run it.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_PY="$ROOT/.venv/bin/python"
say()  { printf '\033[36m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33mWarning: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[31mError: %s\033[0m\n' "$*" >&2; exit 1; }
get_python() {
  if [ -x "$VENV_PY" ]; then echo "$VENV_PY"; return; fi
  for c in python3 python; do command -v "$c" >/dev/null 2>&1 && { command -v "$c"; return; }; done
  return 1
}
set_env() {
  # R packages go to a user library inside the repository, so nothing needs root
  export R_LIBS_USER="$ROOT/R/library"; mkdir -p "$R_LIBS_USER"
  export PYTHONUTF8=1 MPLBACKEND=Agg
  [ -d "$ROOT/tools/linux/bin" ] && case ":$PATH:" in *":$ROOT/tools/linux/bin:"*) ;; *) export PATH="$ROOT/tools/linux/bin:$PATH";; esac
  return 0
}
# "12,13", "1-5", "1..5" or "3 4" -> "3 4" (checked, sorted, unique)
parse_chapters() {
  local spec="${1:-1-13}" out="" part a b
  spec="${spec//,/ }"
  for part in $spec; do
    if [[ $part =~ ^([0-9]+)(-|\.\.)([0-9]+)$ ]]; then a=${BASH_REMATCH[1]}; b=${BASH_REMATCH[3]}; out+=" $(seq "$a" "$b" | tr '\n' ' ')"
    elif [[ $part =~ ^[0-9]+$ ]]; then out+=" $part"
    else die "cannot read the chapter list '$1' (use 12,13 or 1-5)"; fi
  done
  for part in $out; do (( part >= 1 && part <= 13 )) || die "chapters must be between 1 and 13, not $part"; done
  echo $out | tr ' ' '\n' | sort -un | tr '\n' ' ' | sed 's/ *$//'
}
