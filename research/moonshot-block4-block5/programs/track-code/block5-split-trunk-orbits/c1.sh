#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 3000 python3 e1census.py 11 6 120 42 24 3 3 > out_c1_4224_33.txt 2>&1 &
timeout 3000 python3 e1census.py 12 6 120 43 34 4 4 > out_c1_4334_44.txt 2>&1 &
wait
