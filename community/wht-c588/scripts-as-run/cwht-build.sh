#!/bin/bash
# c-wht: overlay the new modules on the copy of Shea's tree and build them, one module per lake call.
# usage: cwht-build.sh <tree root (copy of /content/work/shea)> <input dir with Work/ and comparator/>
set -u
T=$1/lean/complex588; IN=$2; L=/content/work/logs/c-wht
. /content/work/env.sh
cd "$T" || exit 90
echo "=== START $(date -u) tree=$T"
echo "git HEAD: $(git rev-parse HEAD)"
cp -r "$IN"/Work/* Work/
cp "$IN"/comparator/*.json comparator/
sha256sum Work/WHTC588*/*.lean comparator/WHTC588*.json
echo "--- base up-to-date check"
lake build --no-build Work.FourierC588.Main Work.Block.WHT WHTCheck.Solution 2>&1 | tail -3
for m in Work.WHTC588.Challenge Work.WHTC588Neg.Challenge Work.WHTC588.Main Work.WHTC588.Solution; do
  echo "--- build $m  lean procs running: $(pgrep -c lean)"
  s=$(date +%s)
  LEAN_NUM_THREADS=2 /usr/bin/time -v -o $L/build-$m.time lake build $m
  rc=$?
  echo "--- $m rc=$rc wall=$(( $(date +%s)-s ))s  $(grep -E 'Maximum resident' $L/build-$m.time)"
  [ $rc = 0 ] || exit $rc
done
echo "--- no-build check"
lake build --no-build Work.WHTC588.Main Work.WHTC588.Solution Work.WHTC588.Challenge Work.WHTC588Neg.Challenge 2>&1 | tail -3
echo "=== END $(date -u)"
