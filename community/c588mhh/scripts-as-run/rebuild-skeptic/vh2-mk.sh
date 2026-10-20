#!/bin/bash
# v-host-2 (skeptic): compile h-L3-build's sources from MY clean copy (/content/work/vh2/src) into MY build dir,
# plain lean, no Lake; dependencies = m-build tree oleans (read only).  Waits while >= 2 Lean procs > 2 GB run.
SRC=/content/work/vh2/src
B=/content/work/vh2/build
P=/content/work/mb/shea/lean/complex588/.lake/packages
MB=/content/work/mb/shea/lean/complex588/.lake/build/lib/lean
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
LOGD=/content/work/logs/v-host-2
TAG=$1; shift
LP=$B:$MB
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:$P/$p/.lake/build/lib/lean; done
cd $SRC || exit 9
for M in "$@"; do
  while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done
  f=$(echo $M | tr . /).lean
  o=$B/$(echo $M | tr . /)
  mkdir -p $(dirname $o)
  LOG=$LOGD/$TAG-$M.log
  echo "START $M $(date -u +%H:%M:%S)" > $LOG
  LEAN_PATH=$LP LEAN_NUM_THREADS=4 /usr/bin/time -v $LEAN -DautoImplicit=false -R . -o $o.olean -i $o.ilean $f >> $LOG 2>&1
  ec=$?
  echo "EXIT $ec $M $(date -u +%H:%M:%S)" >> $LOG
  w=$(grep -E "Elapsed" $LOG | tail -1 | awk '{print $NF}'); r=$(grep -E "Maximum resident" $LOG | awk '{print $NF}')
  echo "$M EXIT $ec wall=$w maxrss_kb=$r sorry=$(grep -ci sorry $LOG) errors=$(grep -c 'error' $LOG)" >> $LOGD/$TAG.summary
  [ $ec -ne 0 ] && break
done
echo DONE >> $LOGD/$TAG.summary
