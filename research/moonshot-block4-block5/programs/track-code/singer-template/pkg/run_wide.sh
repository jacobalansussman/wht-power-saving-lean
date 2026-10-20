#!/bin/bash
# Colab: wider all-T gate samples on ms-fano prefixes (fano_gate.py), three lanes, single-threaded BLAS
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -f libmp_lin.so ] || gcc -O3 -march=native -shared -fPIC mp.c -o libmp_lin.so
lane() { for a in "$@"; do set -- $a; python3 -B fano_gate.py $1 $2 $3 $4 $5 > out_fg_$1_$2_$3_s$5.txt 2>&1; done; }
lane "FE 4 9 2000 51" "FE 4 10 500 52" "FE 4 11 500 53" "FE 4 12 300 54" "FE 2 5 20000 57" &
lane "FS 4 9 1000 55" "FN 4 9 1000 56" "FE 3 7 10000 58" "FE 3 8 2000 62" &
lane "FS 3 7 3000 59" "FN 3 7 3000 60" "FN 2 5 5000 61" "FS 2 5 5000 63" &
wait
echo DONE
