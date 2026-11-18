#!/usr/bin/env bash
# Check Java, BEAST 2.7.7 (with the MODEL_SELECTION package 1.6.2) and AMRFinderPlus 4.2.7 (with its database, release 2026-10-29.1) and install whatever is missing.
# Usage: scripts/setup-heavy.sh [--chapters 6,9] [--skip-beast] [--skip-amrfinder] [--force] [--dry-run]    Ubuntu or WSL2 Ubuntu, x86-64.
#   --chapters   prepare only what these chapters need (6: BEAST 2 and Java; 9: AMRFinderPlus); default: both     --force   download and install again
#   Needs curl, tar and unzip, and BLAST+ and HMMER for the AMRFinderPlus index; scripts/setup.sh --heavy installs them and then calls this script.
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
CHAPTERS=""; SKIP_BEAST=0; SKIP_AMR=0; FORCE=0; DRY=0
while [ $# -gt 0 ]; do case "$1" in
  --chapters) CHAPTERS="${2:?--chapters needs a list}"; shift 2;; --skip-beast) SKIP_BEAST=1; shift;; --skip-amrfinder) SKIP_AMR=1; shift;;
  --force) FORCE=1; shift;; --dry-run) DRY=1; shift;; -h|--help) sed -n 2,5p "$0"; exit 0;; *) die "unknown option $1";; esac; done
BEAST=1; AMR=1
if [ -n "$CHAPTERS" ]; then BEAST=0; AMR=0; for c in $(parse_chapters "$CHAPTERS"); do [ "$c" = 6 ] && BEAST=1; [ "$c" = 9 ] && AMR=1; done; fi
[ "$SKIP_BEAST" = 0 ] || BEAST=0; [ "$SKIP_AMR" = 0 ] || AMR=0
if [ "$BEAST" = 0 ] && [ "$AMR" = 0 ]; then echo "Nothing to do: Chapters 6 and 9 are the ones that use BEAST 2 and AMRFinderPlus."; exit 0; fi

BEAST_VER=2.7.7; MS_VER=1.6.2; AMR_VER=4.2.7; DB_REL=2026-10-29.1
BEAST_URL="https://github.com/CompEvol/beast2/releases/download/v$BEAST_VER/BEAST.v$BEAST_VER.Linux.x86.tgz"
MS_URL="https://github.com/BEAST2-Dev/model-selection/releases/download/v1.6.0/MODEL_SELECTION.addon.v$MS_VER.zip"
AMR_URL="https://github.com/ncbi/amr/releases/download/amrfinder_v$AMR_VER/amrfinder_binaries_v$AMR_VER.tar.gz"
DB_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/$DB_REL"
DB_FILES="
  AMR.LIB AMRProt-mutation.tsv AMRProt-suppress.tsv
  AMRProt-susceptible.fa AMRProt-susceptible.tsv AMRProt.fa
  AMR_CDS.fa AMR_DNA-Acinetobacter_baumannii.fa AMR_DNA-Acinetobacter_baumannii.tsv
  AMR_DNA-Bordetella_pertussis.fa AMR_DNA-Bordetella_pertussis.tsv AMR_DNA-Campylobacter.fa
  AMR_DNA-Campylobacter.tsv AMR_DNA-Clostridioides_difficile.fa AMR_DNA-Clostridioides_difficile.tsv
  AMR_DNA-Enterococcus_faecalis.fa AMR_DNA-Enterococcus_faecalis.tsv AMR_DNA-Enterococcus_faecium.fa
  AMR_DNA-Enterococcus_faecium.tsv AMR_DNA-Escherichia.fa AMR_DNA-Escherichia.tsv
  AMR_DNA-Klebsiella_oxytoca.fa AMR_DNA-Klebsiella_oxytoca.tsv AMR_DNA-Klebsiella_pneumoniae.fa
  AMR_DNA-Klebsiella_pneumoniae.tsv AMR_DNA-Neisseria_gonorrhoeae.fa AMR_DNA-Neisseria_gonorrhoeae.tsv
  AMR_DNA-Salmonella.fa AMR_DNA-Salmonella.tsv AMR_DNA-Staphylococcus_aureus.fa
  AMR_DNA-Staphylococcus_aureus.tsv AMR_DNA-Streptococcus_pneumoniae.fa AMR_DNA-Streptococcus_pneumoniae.tsv
  ReferenceGeneCatalog.txt ReferenceGeneHierarchy.txt amr_targets.fa
  changelog.txt changes.txt database_format_version.txt
  fam.tsv mapgenelist.txt taxgroup.tsv
  version.txt
"
BIN="$ROOT/tools/linux"; AMRF="$ROOT/amrf"; MS_DIR="$HOME/.beast/2.7/MODEL_SELECTION"
SUDO=""; [ "$(id -u)" = 0 ] || SUDO="sudo "
PROBLEMS=""
ok()   { printf '  \033[32mok\033[0m       %s\n' "$*"; }
need() { printf '  \033[33mmissing\033[0m  %s\n' "$*"; }
fail() { printf '  \033[31mFAILED\033[0m   %s\n' "$*"; PROBLEMS="$PROBLEMS
  - $*"; }
act()  { say "$*"; [ "$DRY" = 1 ] && { echo "    (dry run: not executed)"; return 0; }; return 1; }   # true in a dry run: skip the step
fetch() { # URL DEST: download to a temporary name, then move, so an interrupted download is never mistaken for a finished one
  curl -fL --retry 3 --retry-delay 3 -sS -o "$2.part" "$1" && mv "$2.part" "$2"; }

[ "$(uname -s)" = Linux ] && [ "$(uname -m)" = x86_64 ] || die "the BEAST 2 and AMRFinderPlus builds used here are for x86-64 Linux; this machine is $(uname -s) $(uname -m). Install them by hand (docs/HEAVY-ANALYSES.md)."
for t in curl tar unzip; do command -v "$t" >/dev/null 2>&1 || die "$t is missing: run scripts/setup.sh --heavy (or sudo apt-get install -y curl unzip)"; done
space() { # KB: stop early, with a clear message, when a download would not fit (only asked when something has to be downloaded)
  [ "$DRY" = 1 ] && return 0; local free; free="$(df -Pk "$ROOT" | awk 'NR==2{print $4}')"
  [ "$free" -ge "$1" ] || die "not enough disk space in $ROOT: about $(($1 / 1000)) MB are needed for this download, $((free / 1000)) MB are free"; }

java_major() { local v; v="$(java -version 2>&1 | sed -n '1s/.*version "\([0-9]*\).*/\1/p')"; echo "${v:-0}"; }
mkdir -p "$BIN/bin"

# ---------------------------------------------------------------- BEAST 2 and Java
if [ "$BEAST" = 1 ]; then
  echo "BEAST 2 $BEAST_VER (Chapter 6)"
  beast_cmd() { for c in "$BIN/beast/bin/beast" "$(command -v beast 2>/dev/null || true)"; do
      [ -n "$c" ] && [ -x "$c" ] && "$c" -version 2>&1 | grep -q "v$BEAST_VER" && { echo "$c"; return 0; }; done; return 1; }
  if [ "$FORCE" = 0 ] && B="$(beast_cmd)"; then ok "BEAST 2 $BEAST_VER: $B"
  else
    need "BEAST 2 $BEAST_VER"
    if ! act "download and unpack BEAST 2 $BEAST_VER (95 MB) into tools/linux/beast"; then
      space 360000; rm -rf "$BIN/beast"; fetch "$BEAST_URL" "$BIN/beast.tgz" && tar -xzf "$BIN/beast.tgz" -C "$BIN" && rm -f "$BIN/beast.tgz" || fail "BEAST 2 could not be downloaded or unpacked from $BEAST_URL"
      for t in beast beauti treeannotator packagemanager loganalyser; do [ -f "$BIN/beast/bin/$t" ] && ln -sf "$BIN/beast/bin/$t" "$BIN/bin/$t"; done
      if B="$(beast_cmd)"; then ok "BEAST 2 $BEAST_VER installed: $B"; else fail "BEAST 2 was unpacked but does not start (see: $BIN/beast/bin/beast -version)"; fi
    fi
  fi
  # Java: the BEAST 2 download brings its own; otherwise a system Java 17 or newer is needed
  if [ -x "$BIN/beast/jre/bin/java" ]; then ok "Java: the BEAST 2 download includes its own ($("$BIN/beast/jre/bin/java" -version 2>&1 | head -1))"
  elif command -v java >/dev/null 2>&1 && [ "$(java_major)" -ge 17 ]; then ok "Java $(java_major): $(command -v java)"
  else
    need "Java 17 or newer"
    if ! act "install Java (apt: openjdk-17-jre-headless)"; then
      command -v apt-get >/dev/null 2>&1 || fail "no apt-get here: install Java 17 or newer by hand"
      { ${SUDO}apt-get update -qq && ${SUDO}env DEBIAN_FRONTEND=noninteractive apt-get install -y openjdk-17-jre-headless; } >/dev/null 2>&1 || true
      if command -v java >/dev/null 2>&1 && [ "$(java_major)" -ge 17 ]; then ok "Java $(java_major) installed"; else fail "Java 17 or newer could not be installed (try: sudo apt-get install -y openjdk-17-jre-headless)"; fi
    fi
  fi
  ms_version() { sed -n 's/.*<package[^>]*version=.\([0-9.]*\).*/\1/p' "$MS_DIR/version.xml" 2>/dev/null | head -1; }   # version.xml uses single quotes
  if [ "$FORCE" = 0 ] && [ "$(ms_version)" = "$MS_VER" ]; then ok "MODEL_SELECTION package $MS_VER: $MS_DIR"
  else
    need "MODEL_SELECTION package $MS_VER"
    if ! act "install the MODEL_SELECTION package $MS_VER into ~/.beast/2.7"; then
      mkdir -p "$MS_DIR"; fetch "$MS_URL" "$BIN/ms.zip" && unzip -oq "$BIN/ms.zip" -d "$MS_DIR" && rm -f "$BIN/ms.zip" || fail "the MODEL_SELECTION package could not be downloaded or unpacked from $MS_URL"
      if [ "$(ms_version)" = "$MS_VER" ]; then ok "MODEL_SELECTION package $MS_VER installed"; else fail "MODEL_SELECTION is unpacked in $MS_DIR but its version is '$(ms_version)', not $MS_VER"; fi
    fi
  fi
fi

# ---------------------------------------------------------------- AMRFinderPlus
if [ "$AMR" = 1 ]; then
  echo "AMRFinderPlus $AMR_VER and its database $DB_REL (Chapter 9)"
  amr_ok() { [ -x "$AMRF/amrfinder" ] && [ "$("$AMRF/amrfinder" --version 2>/dev/null | head -1)" = "$AMR_VER" ]; }
  if [ "$FORCE" = 0 ] && amr_ok; then ok "AMRFinderPlus $AMR_VER: $AMRF/amrfinder"
  else
    need "AMRFinderPlus $AMR_VER"
    if ! act "download and unpack the AMRFinderPlus $AMR_VER binaries (6 MB) into amrf/"; then
      space 20000; mkdir -p "$AMRF"; fetch "$AMR_URL" "$AMRF/amrf.tgz" && tar -xzf "$AMRF/amrf.tgz" -C "$AMRF" && rm -f "$AMRF/amrf.tgz" || fail "AMRFinderPlus could not be downloaded or unpacked from $AMR_URL"
      if amr_ok; then ok "AMRFinderPlus $AMR_VER installed: $AMRF/amrfinder"; else fail "AMRFinderPlus was unpacked but does not report version $AMR_VER"; fi
    fi
  fi
  # BLAST+ and HMMER build the database index
  if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v hmmpress >/dev/null 2>&1; then
    need "BLAST+ and HMMER (they build the database index)"
    if ! act "install BLAST+ and HMMER (apt: ncbi-blast+ hmmer)"; then
      { ${SUDO}apt-get update -qq && ${SUDO}env DEBIAN_FRONTEND=noninteractive apt-get install -y ncbi-blast+ hmmer; } >/dev/null 2>&1 || true
      if command -v makeblastdb >/dev/null 2>&1 && command -v hmmpress >/dev/null 2>&1; then ok "BLAST+ and HMMER installed"; else fail "BLAST+ and HMMER are missing (sudo apt-get install -y ncbi-blast+ hmmer)"; fi
    fi
  else ok "BLAST+ and HMMER: $(command -v makeblastdb)"; fi
  db_files_ok() { local f; for f in $DB_FILES; do [ -s "$AMRF/db/$f" ] || return 1; done; }
  db_version() { head -1 "$AMRF/db/version.txt" 2>/dev/null || true; }
  db_ready() { db_files_ok && [ "$(db_version)" = "$DB_REL" ] && [ -f "$AMRF/db/.indexed-$DB_REL" ]; }
  if [ "$FORCE" = 0 ] && db_ready; then ok "AMRFinderPlus database $DB_REL, indexed: $AMRF/db"
  else
    need "AMRFinderPlus database $DB_REL (43 files, 132 MB, indexed)"
    if ! act "download the missing database files from NCBI into amrf/db, then index them"; then
      space 300000; mkdir -p "$AMRF/db"; rm -f "$AMRF/db/.indexed-"*; n=0
      for f in $DB_FILES; do
        if [ "$FORCE" = 1 ] || [ ! -s "$AMRF/db/$f" ]; then fetch "$DB_URL/$f" "$AMRF/db/$f" && n=$((n + 1)) || { fail "could not download $DB_URL/$f"; break; }; fi
      done
      echo "    downloaded $n file(s)"
      if db_files_ok && [ "$(db_version)" = "$DB_REL" ]; then
        if [ -x "$AMRF/amrfinder_index" ] && "$AMRF/amrfinder_index" "$AMRF/db" >/dev/null 2>&1; then touch "$AMRF/db/.indexed-$DB_REL"; ok "AMRFinderPlus database $DB_REL installed and indexed"
        else fail "the database was downloaded but amrf/amrfinder_index failed (are BLAST+ and HMMER installed?)"; fi
      else fail "the database in $AMRF/db is incomplete or is not release $DB_REL (found '$(db_version)')"; fi
    fi
  fi
fi

echo
if [ -n "$PROBLEMS" ]; then warn "some parts could not be installed:$PROBLEMS"; exit 1; fi
if [ "$DRY" = 1 ]; then echo "Dry run: nothing was changed."; else echo "Done."
  [ "$BEAST" = 0 ] || echo "BEAST 2 programs: $BIN/bin (the other scripts add this folder to PATH themselves)"; [ "$AMR" = 0 ] || echo "AMRFinderPlus: $AMRF/amrfinder (database in $AMRF/db)"; fi
