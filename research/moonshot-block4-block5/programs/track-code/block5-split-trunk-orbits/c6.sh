#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 2400 python3 lcensus.py 41 6 300 42 24 > out_c6_4224.txt 2>&1 &
timeout 2400 python3 lcensus.py 42 6 300 42 24 8 > out_c6_4224_dec8.txt 2>&1 &
wait
timeout 2400 python3 lcensus.py 43 6 300 0 0 > out_c6_any.txt 2>&1 &
timeout 2400 python3 lcensus.py 44 6 300 43 34 > out_c6_4334.txt 2>&1 &
wait
