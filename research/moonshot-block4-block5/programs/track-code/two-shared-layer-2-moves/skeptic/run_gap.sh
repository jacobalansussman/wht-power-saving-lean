#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
gcc -O2 -shared -fPIC their/lib/mp.c -o their/lib/libmp_lin.so 2>/dev/null
python3 gapcheck.py cross 6 501 502 503 504 601 602 603 604 > out_gap_crossA.txt 2>&1 &
( python3 gapcheck.py cross 5 31 32 33 34 > out_gap_crossB.txt 2>&1; python3 gapcheck.py sig1 5 41 42 43 44 > out_gap_sig1.txt 2>&1; python3 gapcheck.py sig0 5 61 62 63 > out_gap_sig0.txt 2>&1 ) &
wait
