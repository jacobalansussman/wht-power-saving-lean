#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 gram.py chain3 sort 3 4 5 6 7 > out_chain3.txt 2>&1
python3 analyze.py chain3 5 3 4 5 6 7 > out_analyze_chain3.txt 2>&1
python3 paid.py two 2 3 4 > out_paid_two.txt 2>&1
python3 analyze_paid.py two 3 2 3 4 > out_analyze_paid_two.txt 2>&1
echo DONE >> out_chain3.txt
