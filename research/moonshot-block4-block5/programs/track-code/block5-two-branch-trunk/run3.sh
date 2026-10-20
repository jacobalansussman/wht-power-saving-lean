#!/bin/bash
# frustrated same-frame factors: output 1 only (n1 = 2 array-1 families, n2 = 0), nfac 1 vs control 0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 1 2 3; do timeout 2400 python3 split3.py $s 1 3 54321 1.0,0.7 3000 2 0 > out_r3_s${s}_f1.txt 2>&1 & done
for s in 1 2; do timeout 2400 python3 split3.py $s 0 3 54321 1.0,0.7 3000 2 0 > out_r3_s${s}_f0.txt 2>&1 & done
timeout 2400 python3 split3.py 4 1 4 54321 - 3000 2 2 > out_r3_s4_f1_full.txt 2>&1 &
wait
