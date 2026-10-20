#!/bin/bash
# h-L3-build: like hl3-mk.sh, but pops modules one at a time from a queue file (editable while running).
# usage: hl3-q.sh <tag> <queue file>
T=/content/work/hl3/shea/lean/complex588
P=/content/work/mb/shea/lean/complex588/.lake/packages
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
LOGD=/content/work/logs/h-L3-build
TAG=$1; Q=$2
LP=/dev/shm/hl3/tree
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:/dev/shm/hl3/pk/$p; done
cd $T || exit 9
while true; do
  M=$(head -1 $Q | tr -d ' ')
  [ -z "$M" ] && break
  while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done
  sed -i 1d $Q
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
