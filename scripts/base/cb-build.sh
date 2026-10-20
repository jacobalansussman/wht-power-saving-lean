#!/bin/bash
# c-base: build Shea's lean/complex588 (commit 5239484) with the commands of his scripts/build_and_check.sh,
# each wrapped in /usr/bin/time -v. Departure: LEAN_NUM_THREADS=2 (his default is 1); his script allows the override.
. /content/work/env.sh
P=/content/work/shea/lean/complex588
L=/content/work/logs/cb-build
mkdir -p $L
cd $P || exit 90
export LEAN_NUM_THREADS=2
/content/work/cb-sampler.sh $L/mem-sampler.log & SPID=$!
step() { # step <name> <command...>
  local n=$1; shift
  local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1
  local rc=$?
  local e=$(date +%s)
  echo "=== STEP $n END rc=$rc wall=$((e-s))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"
  return $rc
}
T0=$(date +%s)
echo "BUILD START $(date -u) epoch=$T0 HEAD=$(git rev-parse HEAD) status_lines=$(git status --short | wc -l)"
step cache lake exe cache get
step build lake build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms
BRC=$?
echo "BUILD rc=$BRC"
step axioms lake env lean Work/FourierC588/Axioms.lean
step compare-fourier lake env lean --run tools/Compare.lean Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
step compare-negative lake env lean --run tools/Compare.lean Work.FourierC588Neg.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
echo "(negative control must be nonzero)"
step ours-compare-fourier lake env lean --run /content/work/ours/Compare.lean Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
step ours-compare-negative lake env lean --run /content/work/ours/Compare.lean Work.FourierC588Neg.UniformFourierChallenge \
  Work.FourierC588.Main OAI.PowerSaving.transform_mainC588 OAI.PowerSaving.convolution_mainC588
echo "(negative control must be nonzero)"
step nobuild lake build --no-build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge \
  Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.FourierC588.Axioms
echo "olean count: $(find .lake/build/lib/lean -name '*.olean' | wc -l)  lean sources: $(find OAI Work WHTCheck -name '*.lean' | wc -l)"
echo "git status lines after: $(git status --short | wc -l)"; git status --short | head
du -sh .lake/build .lake/packages
kill $SPID
T1=$(date +%s)
echo "BUILD END $(date -u) total_wall=$((T1-T0))s"
