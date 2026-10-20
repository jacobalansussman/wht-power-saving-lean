#!/bin/bash
# peeling + pollution filter that accepts the trunk's own escapes E2-E4 (Lemma Q)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in 12 14 15 16; do timeout 2400 python3 split.py $s 0 0 3 3 0.7 6 40 54321 30000 3 7 0 pe > out_r9_s${s}_33.txt 2>&1 & done
for s in 32 37 38; do timeout 2400 python3 split.py $s 42 24 3 3 0.7 6 40 54321 30000 3 7 0 pe > out_r9_s${s}_4224.txt 2>&1 & done
for s in 41 45; do timeout 2400 python3 split.py $s 43 34 4 4 0.7 6 40 54321 30000 4 7 0 pe > out_r9_s${s}_4334.txt 2>&1 & done
wait
