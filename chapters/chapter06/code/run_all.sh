#!/bin/bash
# Runs every BEAST 2 analysis of Hands-on 5, one after the other (1 CPU). Idempotent: finished runs are skipped.
# BEAST 2 comes from scripts/setup-heavy.sh (tools/linux/bin) or from your own PATH
export PATH="$(cd "$(dirname "$0")/../../.." && pwd)/tools/linux/bin:$PATH"
command -v beast >/dev/null 2>&1 || { echo "beast not found: run scripts/setup-heavy.sh --chapters 6 from the repository folder" >&2; exit 1; }
cd "$(dirname "$0")/xml" || exit 1
run() { # name xml seed
  if [ -f "done_$1" ]; then echo "skip $1"; return 1; fi
  s=$(date +%s); echo "start $1 $(date +%H:%M:%S)"
  beast -overwrite -seed $3 -threads 1 $2 > "out_$1.txt" 2>&1; echo "end $1 exit $? $(( $(date +%s) - s )) s"; touch "done_$1"; return 0
}
run mcmc_het_strict_s1 mcmc_het_strict.xml 1 && cp mcmc_het_strict.log mcmc_het_strict_s1.log      # copy the log only when the run really ran
run mcmc_het_strict_s2 mcmc_het_strict.xml 2 && cp mcmc_het_strict.log mcmc_het_strict_s2.log
run mcmc_het_relaxed mcmc_het_relaxed.xml 1
for m in het_strict iso dr1 dr2 dr3 het_relaxed; do run ps_$m ps_$m.xml 1; done
echo ALLDONE
