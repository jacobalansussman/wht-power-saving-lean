#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b12.list <<'L'
501 42 24 3 3 e3orb1 0 0 543 7
502 42 24 3 3 e3orb1 0 0 543 7
503 42 24 3 3 e3orb1 0 0 543 7
504 0 0 3 3 e3orb1 0 0 543 7
505 44 44 3 3 e3orb1 0 0 543 7
506 42 24 3 3 e3orb1 0 0 543 8
507 42 24 3 3 e3orb1 0 0 543 8
508 0 0 3 3 e3orb1 0 0 543 8
509 42 34 3 4 e3orb1 0 0 543 7
510 43 34 4 4 e3orb1 0 0 543 7
L
cat b12.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b12_s$0_$5_$9.txt 2>&1'
