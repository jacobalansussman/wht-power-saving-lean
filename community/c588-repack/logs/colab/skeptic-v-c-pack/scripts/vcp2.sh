#!/bin/bash
# v-c-pack (skeptic), second attempt: re-run the decisive checks on c-pack's overlay tree with the WHOLE tree
# bind-mounted read-only in a private mount namespace (so lake cannot delete or change anything), and git told
# to trust the 65534-owned package repos (the first attempt, run as root without this, made lake see an
# "URL change" and delete .lake/packages/mathlib; restored from the base by restore.sh).
# Gated: a heavy step starts only when at most 2 heavy processes (RSS > 3 GB) run and MemAvailable >= need.
. /content/work/env.sh
TREE=/content/work/cpack
P=$TREE/lean/complex588
L=/content/work/logs/v-c-pack
D=/content/work/v-c-pack
mkdir -p $L
export LEAN_NUM_THREADS=1
export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0='*'
RO=/content/work/v-c-pack/ro.sh
heavy() { ps -eo rss=,comm= | awk '$1>3000000 && ($2 ~ /^lean$|nanoda|lean4export/) {n++} END {print n+0}'; }
avail() { awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo; }
gate() { local need=$1; local t0=$(date +%s)
  while :; do h=$(heavy); a=$(avail)
    if [ "$h" -le 2 ] && [ "$a" -ge "$need" ]; then echo "gate ok after $(( $(date +%s)-t0 ))s heavy=$h avail=${a}G"; return 0; fi
    sleep 10; done; }
step() { local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"; }
cd $P || exit 90
echo "VCP2 START $(date -u) HEAD=$(git -c safe.directory='*' rev-parse HEAD) status_lines=$(git status --short | wc -l)"
echo "ro test: $($RO /bin/sh -c 'touch .vcp-rotest 2>&1 || echo read-only-ok')"
echo "mathlib git: $(git -C .lake/packages/mathlib rev-parse HEAD) $(git -C .lake/packages/mathlib remote get-url origin)"
R="Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms"
step nobuild $RO lake build --no-build $R
grep -q "URL has changed" $L/nobuild.log && { echo "ABORT: lake wants to re-clone"; exit 3; }
gate 14; step probe-main $RO lake env lean $D/ProbeMain.lean
gate 10; step probe-chal $RO lake env lean $D/ProbeChal.lean
gate 24; step compare-pos $RO lake env lean --run tools/Compare.lean Work.FourierC588.UniformFourierChallenge Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
gate 24; step compare-neg $RO lake env lean --run tools/Compare.lean Work.FourierC588Neg.UniformFourierChallenge Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
echo "(compare-neg must be nonzero)"
echo "VCP2 END $(date -u) status_lines=$(git status --short | wc -l) oleans=$(find .lake/build/lib/lean -name '*.olean' | wc -l)"
