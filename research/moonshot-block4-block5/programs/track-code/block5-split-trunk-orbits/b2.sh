#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b2.list <<'L'
41 42 24 3 3 esc 4 3 543 7
42 43 34 4 4 esc 4 3 543 7
43 0 0 3 3 esc 4 3 543 7
44 43 24 4 3 esc 4 3 543 7
45 44 44 3 3 esc 4 3 543 7
46 42 24 3 3 esc 4 3 543 8
47 42 34 3 4 esc 4 3 543 7
48 0 0 0 0 esc 4 3 543 7
L
cat b2.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b2_s$0_$5_$9.txt 2>&1'
