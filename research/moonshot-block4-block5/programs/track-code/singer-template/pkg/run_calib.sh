#!/bin/sh
# Colab: build the mod-p library if needed, then the exact end-to-end calibration runs (plain-cost FE prefixes)
cd "$(dirname "$0")"
[ -f libmp_lin.so ] || gcc -O3 -march=native -shared -fPIC mp.c -o libmp_lin.so
python3 -B calib_exact.py 2 6 7 3 > out_calib_2_6_s7.txt 2>&1
python3 -B calib_exact.py 3 9 31 2 > out_calib_3_9_s31.txt 2>&1
echo DONE
