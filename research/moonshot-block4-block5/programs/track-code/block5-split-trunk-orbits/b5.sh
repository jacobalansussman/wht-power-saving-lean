#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b5.list <<'L'
31 42 24 3 3 t2 6 3 543 7 300
32 42 24 3 3 t2 6 3 543 7 300
33 42 24 3 3 t2 6 3 543 7 300
34 0 0 3 3 t2 6 3 543 7 300
35 42 24 3 3 t2 6 3 543 8 300
36 42 24 3 3 t2 6 3 543 8 300
38 44 44 3 3 t2 6 3 543 7 300
39 42 24 3 3 t2 6 3 543 7 300
L
cat b5.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} > out_b5_s$0_$5_$9.txt 2>&1'
