#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b1.list <<'L'
31 42 24 3 3 t2 8 3 543 7
32 42 24 3 3 t2 8 3 543 7
33 42 24 3 3 t2 8 3 543 7
34 0 0 3 3 t2 8 3 543 7
35 42 24 3 3 t2 8 3 543 8
36 42 24 3 3 t2 8 3 543 8
37 43 24 4 3 t2 8 3 543 7
38 44 44 3 3 t2 8 3 543 7
L
cat b1.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b1_s$0_$5_$9.txt 2>&1'
