#!/bin/bash
# designed off-S_Z companions on the close branches, lam (4,4)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 1800 python3 split4.py 51 1 3 0 10 54321 3 > out_r5_s51_c1.txt 2>&1 &
timeout 1800 python3 split4.py 52 2 3 0 10 54321 3 > out_r5_s52_c2.txt 2>&1 &
timeout 1800 python3 split4.py 53 1 4 2 12 54321 3 > out_r5_s53_c1r2.txt 2>&1 &
timeout 1800 python3 split4.py 54 2 4 2 12 54321 3 > out_r5_s54_c2r2.txt 2>&1 &
wait
