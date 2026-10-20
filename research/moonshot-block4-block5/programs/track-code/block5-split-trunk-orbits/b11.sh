#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b11.list <<'L'
401 42 24 3 3 e3orb 0 0 543 7
402 42 24 3 3 e3orb 0 0 543 7
403 42 24 3 3 e3orb 0 0 543 7
404 0 0 3 3 e3orb 0 0 543 7
405 44 44 3 3 e3orb 0 0 543 7
406 42 24 3 3 e3orb 0 0 543 8
407 42 24 3 3 e3orb 0 0 543 8
408 0 0 3 3 e3orb 0 0 543 8
409 42 34 3 4 e3orb 0 0 543 7
410 43 34 4 4 e3orb 0 0 543 7
L
cat b11.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b11_s$0_$5_$9.txt 2>&1'
