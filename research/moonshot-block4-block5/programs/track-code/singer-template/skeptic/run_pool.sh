#!/bin/bash
# v-mw2-singer4: pool-cap test (12-class cap vs full pool vs no one_move filter) on FE prefixes at target cost
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -f libmp_lin.so ] || gcc -O3 -march=native -shared -fPIC mp.c -o libmp_lin.so
python3 -B pooltest.py 2 5 1500 901 > out_pool_2_5.txt 2>&1 &
python3 -B pooltest.py 2 5 1500 902 > out_pool_2_5b.txt 2>&1 &
python3 -B pooltest.py 3 7 800 903 > out_pool_3_7.txt 2>&1 &
python3 -B pooltest.py 4 9 300 904 > out_pool_4_9.txt 2>&1 &
python3 -B pooltest.py 3 8 400 905 > out_pool_3_8.txt 2>&1 &
python3 -B pooltest.py 2 6 300 906 > out_pool_2_6.txt 2>&1 &
wait
echo DONE
