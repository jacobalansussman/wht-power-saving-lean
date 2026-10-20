#!/bin/bash
# h-L3-build: final checks of an instance NAME (C588MHH or C588MH) in MY tree, LEAN_PATH = RAM-disk oleans (hl3-shm.sh):
# (compare part only) Compare.lean (ours, = the tree's tools/Compare.lean) positive and negative, then the official comparator (jail2.sh)
# positive and negative.  usage: hl3-final.sh NAME     (logs: /content/work/logs/h-L3-build/final-NAME-*)
N=$1; T=/content/work/hl3/shea/lean/complex588; L=/content/work/logs/h-L3-build
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
LP=/dev/shm/hl3/tree
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do [ -d /dev/shm/hl3/pk/$p ] && LP=$LP:/dev/shm/hl3/pk/$p; done
cd $T || exit 90
export LEAN_NUM_THREADS=2
export PATH=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin:$PATH
waitslot() { while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done; }
step() { local n=$1; shift; local s=$(date +%s); waitslot
  echo "=== STEP $n START $(date -u +%H:%M:%S)"
  /usr/bin/time -v -o $L/final-$N-$n.time "$@" > $L/final-$N-$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/final-$N-$n.time) last: $(tail -1 $L/final-$N-$n.log)"
  return $rc; }
echo "Compare.lean: $(sha256sum /content/work/ours/Compare.lean tools/Compare.lean | tr '\n' ' ')"
step compare-pos env LEAN_PATH=$LP $LEAN --run /content/work/ours/Compare.lean Work.Fourier$N.UniformFourierChallenge \
  Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N
step compare-neg env LEAN_PATH=$LP $LEAN --run /content/work/ours/Compare.lean Work.Fourier${N}Neg.UniformFourierChallenge \
  Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N
echo "(negative comparison must be nonzero: rc=$?)"
grep -E "RESULT" $L/final-$N-compare-pos.log; grep -E "definition body differs|RESULT" $L/final-$N-compare-neg.log | head -3
echo COMPARE DONE
