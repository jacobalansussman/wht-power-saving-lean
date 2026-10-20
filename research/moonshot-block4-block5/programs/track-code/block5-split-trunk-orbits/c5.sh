#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 3000 python3 cover.py 31 6 400 42 24 3 3 > out_c5_4224_33.txt 2>&1 &
timeout 3000 python3 cover.py 32 8 300 0 0 0 0 > out_c5_any.txt 2>&1 &
wait
