#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b10.list <<'L'
301 42 24 3 3 side 0 1 543 7
302 42 24 3 3 side 0 1 543 7
303 0 0 3 3 side 0 1 543 7
304 44 44 3 3 side 0 1 543 7
305 42 24 3 3 sideo 0 2 543 7
306 42 24 3 3 side 0 1 543 8
307 42 24 3 3 sideo 0 2 543 8
308 0 0 0 0 sideo 0 2 543 7
L
cat b10.list | xargs -P 2 -L 1 bash -c 'timeout 600 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b10_s$0_$5_$9.txt 2>&1'
