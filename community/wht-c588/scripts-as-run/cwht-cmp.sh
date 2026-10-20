#!/bin/bash
# c-wht: official comparator runs through c-base's cb-cmp.sh (jail2.sh stand-in sandbox), one at a time.
# usage: cwht-cmp.sh <tree root> <label prefix> <config>...    (configs relative to lean/complex588)
set -u
P=$1/lean/complex588; PFX=$2; shift 2
L=/content/work/logs/c-wht/cmp; mkdir -p $L
waitslot() { local w=0; while [ $(pgrep -c -x lean) -ge 2 ]; do sleep 10; w=$((w+10)); done; echo "  (slot free after waiting ${w}s; lean procs now: $(pgrep -c -x lean))"; }
echo "=== chown .lake to 65534 $(date -u)"; s=$(date +%s)
chown -R 65534:65534 $P/.lake; echo "chown rc=$? $(( $(date +%s)-s ))s"
for c in "$@"; do
  lab=$PFX-$(basename $c .json)
  echo "=== $lab $(date -u)"; waitslot
  /content/work/cb-cmp.sh $lab $P $c > $L/$lab.log 2>&1; rc=$?
  cp /content/work/logs/cb-cmp/$lab.time $L/ 2>/dev/null
  echo "rc=$rc"; grep -E "accepts|okay|does not match|Const|error|Error|END rc|Maximum resident|Elapsed" $L/$lab.log | cut -c1-300
done
echo "=== END $(date -u)"
