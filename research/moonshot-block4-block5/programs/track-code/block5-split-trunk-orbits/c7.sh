#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 2400 python3 orbstat.py 61 6 40000 42 24 3 3 > out_c7_4224.txt 2>&1 &
timeout 2400 python3 orbstat.py 62 6 40000 42 24 3 3 8 > out_c7_4224_dec8.txt 2>&1 &
wait
timeout 2400 python3 orbstat.py 63 6 40000 0 0 0 0 > out_c7_any.txt 2>&1 &
timeout 2400 python3 orbstat.py 64 6 40000 0 0 0 0 8 > out_c7_any_dec8.txt 2>&1 &
wait
