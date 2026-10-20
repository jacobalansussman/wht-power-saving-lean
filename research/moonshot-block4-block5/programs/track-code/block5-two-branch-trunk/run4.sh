#!/bin/bash
# lam (4,4) on the CLOSE branches (profiles (4,3) / (3,4)): the only lam=4 door left by Lemma L4'; + census robustness
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 41 42 43 44; do timeout 2400 python3 split.py $s 43 34 4 4 0.7 3 30 54321 4000 4 7 0 > out_r4_s${s}.txt 2>&1 & done
timeout 2400 python3 cfgprobe2.py 7 30 > out_cfgprobe2_s7_n30.txt 2>&1 &
wait
