#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 paid.py cliff2 2 3 4 > out_paid_cliff2.txt 2>&1
python3 analyze_paid.py cliff2 3 2 3 4 > out_analyze_paid_cliff2.txt 2>&1
python3 paid.py h1small 2 3 4 5 > out_paid_h1small.txt 2>&1
python3 analyze_paid.py h1small 3 2 3 4 5 > out_analyze_paid_h1small.txt 2>&1
echo DONE >> out_paid_h1small.txt
