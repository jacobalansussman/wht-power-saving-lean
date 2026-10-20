#!/bin/bash
# batch 3 on the Mac (one heavy slot): every line of jobs_b3.txt, menus 2 and 3, 12 skeletons x 2 starts
cd "$(dirname "$0")"
PY=python3
mkdir -p b3
while read m W H shape; do
  for mn in 2 3; do
    tag=$(echo "${m}_${W}_${H}_${shape}_mn${mn}" | tr ':,' '-_')
    $PY symsolve.py fast $m $W $H $shape $mn 2 $RANDOM 12 6000 > b3/b3_$tag.txt 2>&1
  done
done < jobs_b3.txt
echo DONE > b3/DONE
