#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in $1; do for seed in 1 2 3; do timeout 1500 python3 loc.py 4 $s $seed 1500 > out_lm_n4_s${s}_seed$seed.txt 2>&1 & done; done; wait
