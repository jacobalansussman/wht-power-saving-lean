#!/bin/bash
# k-lower: compile the new modules ONE AT A TIME with `lake env lean -o` (no Lake graph replay per module: on the loaded
# Colab each `lake build <m>` spent minutes replaying ~9100 traces).  Same options as the lakefile (autoImplicit false),
# LEAN_NUM_THREADS=2; waits while 2 or more HEAVY Lean processes (RSS > 2 GB) run (the 2-process rule; light Lean jobs of other builds, e.g. Aesop files, do not block).  Oleans go to the package's .lake/build/lib/lean, so
# `lake env` sees them; Lake traces are NOT written (a later `lake build` would recompile these modules).
# usage: build-direct.sh <PKG dir> <tree lean/complex588> <log dir> [first module index]
set -uo pipefail
. /content/work/env.sh
PKG=$1; T=$2; L=$3; FIRST=${4:-1}
mkdir -p $L
cd $T
echo "DIRECT BUILD START $(date -u) pkg=$PKG"
n=0
while read m; do
  [ -z "$m" ] && continue
  n=$((n+1))
  [ $n -lt $FIRST ] && continue
  while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done
  f=$(echo $m | tr . /).lean
  o=.lake/build/lib/lean/$(echo $m | tr . /)
  mkdir -p $(dirname $o)
  log=$L/$(printf %02d $n)-$m.log
  LEAN_NUM_THREADS=2 /usr/bin/time -v lake env lean -DautoImplicit=false -R . -o $o.olean -i $o.ilean $f > $log 2>&1
  rc=$?
  t=$(grep -E "Elapsed" $log | tail -1 | awk '{print $NF}'); r=$(grep -E "Maximum resident" $log | awk '{print $NF}')
  ax=$(grep -c "depends on axioms" $log); so=$(grep -c "sorryAx" $log)
  echo "$(printf %02d $n) rc=$rc $m wall=$t maxrss_kb=$r axioms_lines=$ax sorryAx=$so"
  if [ $rc -ne 0 ] || [ $so -ne 0 ]; then echo "FAILED at $m (log $log)"; grep -E "error" $log | head -20; exit 1; fi
done < $PKG/MODULES.txt
echo "DIRECT BUILD END $(date -u)"
