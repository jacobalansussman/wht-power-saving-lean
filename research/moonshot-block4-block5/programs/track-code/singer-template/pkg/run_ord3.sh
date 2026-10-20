#!/bin/bash
# Colab: Singer templates with ORDER-3 filters (v-ms-fano's sigma, S4_ORDER=3) through the all-T gate (t_full.py); 2 lanes
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 S4_ORDER=3
[ -f libmp_lin.so ] || gcc -O3 -march=native -shared -fPIC mp.c -o libmp_lin.so
python3 -B chk_sigma.py > out_o3_chk_sigma.txt 2>&1
lane() { for a in "$@"; do set -- $a; python3 -B t_full.py $1 $2 $3 $4 $5 > out_o3_$5_$1_$2_s$4.txt 2>&1; done; }
lane "2 8 30 80 PL" "2 7 300 81 HWs" "3 10 200 82 HWs" &
lane "2 7 300 83 HWc" "2 7 60 84 HW" "2 7 60 85 H" &
wait
echo DONE
