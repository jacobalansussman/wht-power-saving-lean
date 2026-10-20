#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 2000 python3 e1census.py 21 6 160 42 24 3 3 8 > out_c4_dec8_4224_33.txt 2>&1 &
timeout 2000 python3 e1census.py 22 6 160 0 0 0 0 8 > out_c4_dec8_any.txt 2>&1 &
wait
