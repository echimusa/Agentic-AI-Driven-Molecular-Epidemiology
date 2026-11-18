#!/usr/bin/env bash
# Copy this repository from a Windows drive (/mnt/c/...) into the Linux home folder of WSL2, then run scripts/setup.sh there with the same options.
# Started by scripts/setup.ps1 on Windows, or by scripts/setup.sh when it finds itself on a Windows drive. Run it inside WSL2 Ubuntu.
# Usage: bash scripts/wsl-handoff.sh [options of setup.sh, for example --chapters 5,6,7 --dry-run]
# MEPI_WSL_DIR sets the destination (default ~/<name of the repository folder>); MEPI_SRC names the source folder when this file is run from elsewhere.
set -euo pipefail
SRC="${MEPI_SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DST="${MEPI_WSL_DIR:-$HOME/$(basename "$SRC")}"
say() { printf '\033[36m==> %s\033[0m\n' "$*"; }
DRY=0; for a in "$@"; do [ "$a" = "--dry-run" ] && DRY=1; done
if [ "$DRY" = 1 ]; then
  say "dry run: would copy $SRC to $DST (leaving out .venv, R/library, .work and tools/win) and run: bash scripts/setup.sh --here $*"
  exit 0
fi
if [ "$SRC" != "$DST" ]; then
  say "copy the repository to $DST (leaving out .venv, R/library, .work and tools/win, which are built for each system)"
  mkdir -p "$DST"
  (cd "$SRC" && tar --exclude='.venv' --exclude='.work' --exclude='R/library' --exclude='tools/win' --exclude='__pycache__' -cf - .) | tar -xf - -C "$DST"
  # files saved on Windows can carry CRLF line endings, which bash cannot run
  for f in "$DST"/scripts/*.sh "$DST"/tools/*.py "$DST"/r-packages.R "$DST"/requirements.txt; do [ -f "$f" ] && sed -i 's/\r$//' "$f"; done
fi
# the book was validated on Ubuntu 24.04 with Python 3.12; "wsl --install -d Ubuntu" can give a newer release
. /etc/os-release 2>/dev/null || true
if [ "${VERSION_ID:-}" != "24.04" ]; then
  printf '\033[33mWarning: this is %s with %s; the book was validated on Ubuntu 24.04 with Python 3.12.3. If a package fails to install, see docs/TROUBLESHOOTING.md; "wsl --list --online" shows the Ubuntu releases this PC can install.\033[0m\n' "${PRETTY_NAME:-an unknown Linux}" "$(python3 --version 2>&1 || echo 'no Python yet')" >&2
fi
say "run setup in $DST"
rc=0; (cd "$DST" && bash scripts/setup.sh --here "$@") || rc=$?
echo
say "the Linux copy is $DST"
echo "  Run a chapter:          cd $DST && bash scripts/run-chapter.sh 5"
printf '  Open it from Windows:   \\\\wsl$\\%s%s\n' "${WSL_DISTRO_NAME:-Ubuntu}" "$(printf '%s' "$DST" | tr / '\\')"
echo "  Changes made to the Windows copy reach the Linux copy when you run setup.ps1 (or this script) again."
exit "$rc"
