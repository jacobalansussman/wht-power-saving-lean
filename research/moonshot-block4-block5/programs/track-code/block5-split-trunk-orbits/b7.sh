#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b7.list <<'L'
81 42 24 3 3 t2 0 0 543 7 0 s2b
82 42 24 3 3 t2 0 0 543 7 0 s2b
83 0 0 3 3 t2 0 0 543 7 0 s2b
84 44 44 3 3 t2 0 0 543 7 0 s2b
85 42 24 3 3 t2 0 0 543 7 0 s2b,fixNb
86 0 0 3 3 t2 0 0 543 7 0 s2b,fixNb
87 42 24 3 3 orbit 0 0 543 7 0 s2b,fixNb
88 0 0 3 3 orbit 0 0 543 7 0 s2b,fixNb
89 42 24 3 3 t2 4 3 543 7 300 s2b,fixNb
90 42 24 3 3 t2 0 0 543 8 0 s2b,fixNb
L
cat b7.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} > out_b7_s$0_$5_$9.txt 2>&1'
