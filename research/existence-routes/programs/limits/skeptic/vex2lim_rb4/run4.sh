#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 paid.py cliff2s 2 3 4 5 > out_paid_cliff2s.txt 2>&1
python3 analyze_paid.py cliff2s 4 2 3 4 5 > out_analyze_paid_cliff2s.txt 2>&1
echo DONE >> out_paid_cliff2s.txt
