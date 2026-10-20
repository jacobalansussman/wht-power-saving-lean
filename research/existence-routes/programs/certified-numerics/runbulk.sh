#!/bin/bash
# runbulk.sh: 3 workers over all a2-nag exported relaxations (smallest files first), outputs out/bulk_<k>.jsonl
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ls -S -r a2/*/*.json > out/bulk_files.txt
for k in 0 1 2; do
  awk -v k=$k 'NR%3==k' out/bulk_files.txt > out/bulk_list_$k.txt
  ( xargs -n 20 timeout 3000 python3 bulk.py out/bulk_$k.jsonl < out/bulk_list_$k.txt > out/bulk_$k.log 2>&1; echo DONE >> out/bulk_$k.log ) &
done
wait
