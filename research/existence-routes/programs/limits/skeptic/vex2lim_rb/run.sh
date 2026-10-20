#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 ctl.py > out_ctl.txt 2>&1
python3 gram.py main sort 2 3 4 5 6 > out_main_sort_2to6.txt 2>&1
python3 gram.py main dense 2 3 4 5 6 7 > out_main_dense_2to7.txt 2>&1
python3 analyze.py main 5 3 4 5 6 7 > out_analyze_main.txt 2>&1
python3 paid.py one 2 3 4 5 > out_paid_one.txt 2>&1
python3 analyze_paid.py one 4 3 4 5 > out_analyze_paid_one.txt 2>&1
echo DONE >> out_ctl.txt
