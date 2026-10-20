#!/bin/bash
# split trunk with injection rounds; lam (3,3) and (4,4); other-branch images at layer <= 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 11 12 13; do timeout 1500 python3 split.py $s 0 0 3 3 0.7 4 40 54321 2500 3 > out_r1_s${s}_33.txt 2>&1 & done
for s in 21 22 23; do timeout 1500 python3 split.py $s 0 0 4 4 0.7 4 40 54321 2500 3 > out_r1_s${s}_44.txt 2>&1 & done
wait
