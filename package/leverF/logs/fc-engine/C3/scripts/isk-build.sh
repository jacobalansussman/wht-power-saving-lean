#!/bin/bash
# i-skel: lake build of the given targets in my tree, waiting for a heavy slot first.
# usage: isk-build.sh <TAG> <target> [target ...]      (LEAN_NUM_THREADS from ISK_THREADS, default 3 minus the other heavy Lean processes)
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
TAG=$1; shift
L=$LOGROOT/$(lower $TAG); mkdir -p $L
cd $T || fail "no $T"
[ "$(stat -c %U .lake)" = root ] || { echo "chown .lake back to root"; chown -R root:root .lake; }
wl 12000
# Lake runs about LEAN_NUM_THREADS Lean jobs at once: keep my jobs plus the machine's other heavy Lean processes <= 3
THR=${ISK_THREADS:-$(( MAXHEAVY - $(heavy) ))}; [ $THR -ge 1 ] || THR=1
echo "   LEAN_NUM_THREADS=$THR"
step build env LEAN_NUM_THREADS=$THR lake build "$@"; rc=$?
grep -E "✖|error:|Build completed" $L/build.log | head -30
grep -E "Built .* \([0-9.]+(s|m)\)" $L/build.log | sed -E 's/.*Built ([^ ]+) \(([0-9.]+)(s|m)\).*/   \2\3 \1/' | sort -g -r | head -8
[ $rc -eq 0 ] || fail "build $TAG rc=$rc"
step nobuild lake build --no-build "$@"; tail -1 $L/nobuild.log
echo "BUILD $TAG OK"
