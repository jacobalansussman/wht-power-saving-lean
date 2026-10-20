#!/bin/bash
# run_top.sh: four top-layer diagnostics in parallel (12 GB virtual cap each)
ulimit -v 12000000
for k in plain pure inject sepT; do
  python3 -B topdiag.py $k 6 7 > out_top_$k.txt 2>&1 &
done
wait
