#!/bin/bash
# c-pack: rebuild the overlay tree (only changed modules and their dependents recompile), then the checks.
# Mirrors c-base's cb-build.sh (Shea's scripts/build_and_check.sh commands, each in /usr/bin/time -v).
. /content/work/env.sh
P=/content/work/cpack/lean/complex588
L=/content/work/logs/c-pack/build
mkdir -p $L
cd $P || exit 90
export LEAN_NUM_THREADS=${CP_THREADS:-1}
/content/work/cb-sampler.sh $L/mem-sampler.log & SPID=$!
step() { local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"
  return $rc; }
T0=$(date +%s)
echo "BUILD START $(date -u) LEAN_NUM_THREADS=$LEAN_NUM_THREADS HEAD=$(git rev-parse HEAD) status_lines=$(git status --short | wc -l)"
step build lake build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms
BRC=$?; echo "BUILD rc=$BRC"
[ $BRC = 0 ] || { tail -60 $L/build.log; kill $SPID; exit 1; }
echo "modules rebuilt (Built lines): $(grep -c 'Built ' $L/build.log)   replayed: $(grep -c 'Replayed ' $L/build.log)"
step axioms lake env lean Work/FourierC588/Axioms.lean
step ours-compare-fourier lake env lean --run /content/work/ours/Compare.lean Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
step ours-compare-negative lake env lean --run /content/work/ours/Compare.lean Work.FourierC588Neg.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
echo "(negative control must be nonzero)"
step nobuild lake build --no-build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms
echo "olean count: $(find .lake/build/lib/lean -name '*.olean' | wc -l)  lean sources: $(find OAI Work WHTCheck -name '*.lean' | wc -l)"
kill $SPID
echo "BUILD END $(date -u) total_wall=$(( $(date +%s)-T0 ))s"
