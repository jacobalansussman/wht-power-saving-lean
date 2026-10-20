#!/bin/bash
# peeling closure: controls (t0-type, no injections) + greedy rounds on run1/run2 configurations with exact forced zeros
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 600 python3 split.py 3 0 0 0 0 0.7 0 0 543 2500 4 7 0 p > out_r6_ctrl_s3.txt 2>&1 &
timeout 2400 python3 split.py 12 0 0 3 3 0.7 3 30 54321 20000 3 7 0 p > out_r6_s12_33.txt 2>&1 &
timeout 2400 python3 split.py 32 42 24 3 3 0.7 3 30 54321 20000 3 7 0 p > out_r6_s32_4224.txt 2>&1 &
timeout 2400 python3 split.py 41 43 34 4 4 0.7 3 30 54321 20000 4 7 0 p > out_r6_s41_4334.txt 2>&1 &
timeout 2400 python3 split4.py 51 1 3 2 10 54321 3 0.7 p 20000 > out_r6_s51_comp.txt 2>&1 &
wait
