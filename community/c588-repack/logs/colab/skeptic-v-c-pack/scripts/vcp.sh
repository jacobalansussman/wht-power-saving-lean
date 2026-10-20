#!/bin/bash
# v-c-pack (skeptic): re-run the decisive checks on c-pack's overlay tree. Gated: a heavy step starts only when
# at most 2 other heavy processes (RSS > 3 GB: lean, nanoda, lean4export) run and MemAvailable >= need.
. /content/work/env.sh
P=/content/work/cpack/lean/complex588
L=/content/work/logs/v-c-pack
D=/content/work/v-c-pack
mkdir -p $L
cd $P || exit 90
export LEAN_NUM_THREADS=1
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
echo "VCP START $(date -u) HEAD=$(git rev-parse HEAD) status_lines=$(git status --short | wc -l)"
sha256sum Work/FourierC588/UniformFourierChallenge.lean Work/FourierC588Neg/UniformFourierChallenge.lean Work/FourierC588/Goal.lean Work/FourierC588/Main.lean tools/Compare.lean
R="Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms"
step nobuild lake build --no-build $R
gate 14; step probe-main lake env lean $D/ProbeMain.lean
gate 10; step probe-chal lake env lean $D/ProbeChal.lean
gate 24; step compare-pos lake env lean --run tools/Compare.lean Work.FourierC588.UniformFourierChallenge Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
gate 24; step compare-neg lake env lean --run tools/Compare.lean Work.FourierC588Neg.UniformFourierChallenge Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
echo "(compare-neg must be nonzero)"
step nobuild2 lake build --no-build $R
echo "VCP END $(date -u) status_lines=$(git status --short | wc -l) oleans=$(find .lake/build/lib/lean -name '*.olean' | wc -l)"
