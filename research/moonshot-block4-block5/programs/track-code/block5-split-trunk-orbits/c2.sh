#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 3000 python3 e1census.py 13 6 120 43 24 4 3 > out_c2_4324_43.txt 2>&1 &
timeout 3000 python3 e1census.py 14 6 120 44 44 3 3 > out_c2_4444_33.txt 2>&1 &
wait
