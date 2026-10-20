#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
gcc -O2 -shared -fPIC their/lib/mp.c -o their/lib/libmp_lin.so 2>/dev/null
( python3 pools.py cross 5 31 33 > out_pools_cross.txt 2>&1 ) &
( python3 pools.py sig1 5 43 44 41 42 > out_pools_sig1.txt 2>&1 ) &
wait
