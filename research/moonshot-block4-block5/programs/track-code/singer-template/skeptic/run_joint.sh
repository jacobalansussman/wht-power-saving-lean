#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -f libmp_lin.so ] || gcc -O3 -march=native -shared -fPIC mp.c -o libmp_lin.so
python3 -B jointcheck.py 2 5 400 911 > out_joint_2_5.txt 2>&1 &
python3 -B jointcheck.py 2 5 400 912 > out_joint_2_5b.txt 2>&1 &
python3 -B jointcheck.py 2 6 150 913 > out_joint_2_6.txt 2>&1 &
python3 -B jointcheck.py 3 7 100 914 > out_joint_3_7.txt 2>&1 &
wait
echo DONE
