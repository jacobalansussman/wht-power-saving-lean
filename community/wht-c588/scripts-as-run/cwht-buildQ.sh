#!/bin/bash
# c-wht: overlay the new modules on the copy of Shea's tree and build them, one module per lake call.
# usage: cwht-build.sh <tree root (copy of /content/work/shea)> <input dir with Work/ and comparator/>
set -u
T=$1/lean/complex588; IN=$2; L=/content/work/logs/c-wht/packed; mkdir -p $L
. /content/work/env.sh
cd "$T" || exit 90
waitslot() { local w=0; while [ $(pgrep -c -x lean) -ge 2 ]; do sleep 10; w=$((w+10)); done; echo "  (slot free after waiting ${w}s)"; }
echo "=== START $(date -u) tree=$T"
echo "git HEAD: $(git rev-parse HEAD)"
cp -r "$IN"/Work/* Work/
cp "$IN"/comparator/*.json comparator/
sha256sum Work/WHTC588*/*.lean comparator/WHTC588*.json
echo "--- base up-to-date check (GATE: must print All targets up-to-date, else stop before building anything)"
NB=$(lake build --no-build Work.FourierC588.Main Work.FourierC588.UniformFourierChallenge Work.FourierC588Neg.UniformFourierChallenge Work.GCert.Data.C588Controls Work.Block.WHT WHTCheck.Solution 2>&1)
echo "$NB" | grep -v "^info:" | tail -5
echo "$NB" | grep -q "All targets up-to-date" || { echo "GATE FAILED: tree not up to date, nothing built"; exit 93; }
for m in Work.WHTC588Q.Challenge Work.WHTC588QNeg.Challenge Work.WHTC588Q.Main Work.WHTC588Q.Solution; do
  echo "--- build $m  lean procs running: $(pgrep -c lean)"
  waitslot
  s=$(date +%s)
  LEAN_NUM_THREADS=2 /usr/bin/time -v -o $L/build-$m.time lake build $m
  rc=$?
  echo "--- $m rc=$rc wall=$(( $(date +%s)-s ))s  $(grep -E 'Maximum resident' $L/build-$m.time)"
  [ $rc = 0 ] || exit $rc
done
echo "--- no-build check"
lake build --no-build Work.WHTC588Q.Main Work.WHTC588Q.Solution Work.WHTC588Q.Challenge Work.WHTC588QNeg.Challenge 2>&1 | tail -3
echo "=== END $(date -u)"
