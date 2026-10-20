#!/bin/bash
# split trunk, most separated profiles: branch a (4,2) for array 1, branch b (2,4) for array 2; lam (3,3); LIMY 3; no double
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 31 32 33 34; do timeout 2400 python3 split.py $s 42 24 3 3 0.7 3 40 54321 5000 3 7 0 > out_r2_s${s}.txt 2>&1 & done
for s in 35 36; do timeout 2400 python3 split.py $s 42 24 3 3 0.7 3 40 54321 5000 3 8 0 > out_r2_s${s}_dec8.txt 2>&1 & done
wait
