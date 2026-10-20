#!/bin/bash
# rfwht-cmp.sh <tree root> <config>...   (rf-wht, 13 Oct 2026; r-wht's rwht-cmp.sh with the rf-wht paths, and the memory
# guard of rfwht-run.sh in the place of the 2-Lean-processes wait)
# Official comparator runs, one at a time, through rfwht-cb-cmp.sh (c-base's cb-cmp.sh: jail2.sh stand-in for landrun).
# Builds the private seccomp helper first. After the chown below, do NOT run `lake` as root in this tree.
set -u
P=$1/lean/complex588; shift
L=/content/work/logs/rf-wht/cmp; mkdir -p $L
D=/content/work/rf-wht/cc
gcc -O2 -Wall -o $D/nounix $D/nounix.c && echo "nounix built: $(sha256sum $D/nounix)"
chmod 755 $D $D/nounix $D/jail2.sh $D/rfwht-cb-cmp.sh
memok() { [ $(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo) -ge 16 ]; }
waitmem() { local w=0; while ! memok; do sleep 10; w=$((w+10)); done; echo "  (waited ${w}s for 16 GB MemAvailable; now $(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo) GB)"; }
echo "=== chown .lake to 65534 $(date -u) (not following the packages symlink)"; s=$(date +%s)
chown -R 65534:65534 $P/.lake; echo "chown rc=$? $(( $(date +%s)-s ))s; packages: $(ls -ld $P/.lake/packages)"
for c in "$@"; do
  lab=rfwht-$(basename $c .json)
  echo "=== $lab $(date -u)"; waitmem
  $D/rfwht-cb-cmp.sh $lab $P $c > $L/$lab.log 2>&1; rc=$?
  echo "rc=$rc"; grep -E "accepts|okay|does not match|Const|error|Error|END rc|Maximum resident|Elapsed" $L/$lab.log | cut -c1-300
done
echo "=== END $(date -u)"
