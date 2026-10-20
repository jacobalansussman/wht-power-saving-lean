#!/bin/bash
# h-L3-build: compile modules of MY copy of the tree with plain `lean` (no Lake at all: `lake build` in a copied tree
# re-resolves the Mathlib URL and deletes .lake/packages/mathlib).  LEAN_PATH = my tree's build dir + the package build
# dirs of the m-build tree (read only).  One module at a time, LEAN_NUM_THREADS=2; waits while 2 or more Lean
# processes with RSS > 2 GB run (the 2-heavy-processes rule).
# usage: hl3-mk.sh <tag> <module> [<module> ...]      logs: /content/work/logs/h-L3-build/<tag>-<module>.log
T=/content/work/hl3/shea/lean/complex588
P=/content/work/mb/shea/lean/complex588/.lake/packages
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
LOGD=/content/work/logs/h-L3-build
TAG=$1; shift
LP=$T/.lake/build/lib/lean
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:$P/$p/.lake/build/lib/lean; done
cd $T || exit 9
for M in "$@"; do
  while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done
  f=$(echo $M | tr . /).lean
  o=.lake/build/lib/lean/$(echo $M | tr . /)
  mkdir -p $(dirname $o)
  LOG=$LOGD/$TAG-$M.log
  echo "START $M $(date -u +%H:%M:%S)" > $LOG
  LEAN_PATH=$LP LEAN_NUM_THREADS=2 /usr/bin/time -v $LEAN -DautoImplicit=false -R . -o $o.olean -i $o.ilean $f >> $LOG 2>&1
  ec=$?
  echo "EXIT $ec $M $(date -u +%H:%M:%S)" >> $LOG
  w=$(grep -E "Elapsed" $LOG | tail -1 | awk '{print $NF}'); r=$(grep -E "Maximum resident" $LOG | awk '{print $NF}')
  echo "$M EXIT $ec wall=$w maxrss_kb=$r sorryAx=$(grep -c sorryAx $LOG) errors=$(grep -c 'error' $LOG)" >> $LOGD/$TAG.summary
  [ $ec -ne 0 ] && break
done
echo DONE >> $LOGD/$TAG.summary
