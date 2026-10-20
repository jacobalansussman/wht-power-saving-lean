#!/bin/bash
# greedy exact rounds WITH peeling (targets = kill-trigger rows): exact verdict per round
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 12 14 15; do timeout 3000 python3 split.py $s 0 0 3 3 0.7 6 40 54321 30000 3 7 0 p > out_r7_s${s}_33.txt 2>&1 & done
for s in 32 37; do timeout 3000 python3 split.py $s 42 24 3 3 0.7 6 40 54321 30000 3 7 0 p > out_r7_s${s}_4224.txt 2>&1 & done
for s in 41 45; do timeout 3000 python3 split.py $s 43 34 4 4 0.7 6 40 54321 30000 4 7 0 p > out_r7_s${s}_4334.txt 2>&1 & done
timeout 3000 python3 split4.py 51 1 3 6 30 54321 3 0.7 p 30000 > out_r7_s51_comp.txt 2>&1 &
timeout 3000 python3 split.py 36 42 24 3 3 0.7 6 40 54321 30000 3 8 0 p > out_r7_s36_dec8.txt 2>&1 &
wait
