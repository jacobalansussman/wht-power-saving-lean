#!/bin/bash
# rwht-cmp.sh <tree root> <config>...   (r-wht, 13 Oct 2026; same procedure as c-wht's cwht-cmp.sh)
# Official comparator runs, one at a time, through rwht-cb-cmp.sh (c-base's cb-cmp.sh: jail2.sh stand-in for landrun).
# Builds the private seccomp helper first. After the chown below, do NOT run `lake` as root in this tree
# (m-build's gotcha: Lake then deletes and re-clones .lake/packages).
set -u
P=$1/lean/complex588; shift
L=/content/work/logs/r-wht/cmp; mkdir -p $L
D=/content/work/r-wht/cc
gcc -O2 -Wall -o $D/nounix $D/nounix.c && echo "nounix built: $(sha256sum $D/nounix)"
chmod 755 $D $D/nounix $D/jail2.sh $D/rwht-cb-cmp.sh
waitslot() { local w=0; while [ $(pgrep -c -x lean) -ge 2 ]; do sleep 10; w=$((w+10)); done; echo "  (slot free after waiting ${w}s; lean procs now: $(pgrep -c -x lean))"; }
echo "=== chown .lake to 65534 $(date -u)"; s=$(date +%s)
chown -R 65534:65534 $P/.lake; echo "chown rc=$? $(( $(date +%s)-s ))s"
for c in "$@"; do
  lab=rwht-$(basename $c .json)
  echo "=== $lab $(date -u)"; waitslot
  $D/rwht-cb-cmp.sh $lab $P $c > $L/$lab.log 2>&1; rc=$?
  echo "rc=$rc"; grep -E "accepts|okay|does not match|Const|error|Error|END rc|Maximum resident|Elapsed" $L/$lab.log | cut -c1-300
done
echo "=== END $(date -u)"
