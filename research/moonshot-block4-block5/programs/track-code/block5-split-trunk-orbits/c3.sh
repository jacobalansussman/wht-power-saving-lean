#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 3000 python3 e1census.py 15 6 120 0 0 0 0 > out_c3_any.txt 2>&1 &
timeout 3000 python3 e1census.py 16 6 120 42 34 3 4 > out_c3_4234_34.txt 2>&1 &
wait
