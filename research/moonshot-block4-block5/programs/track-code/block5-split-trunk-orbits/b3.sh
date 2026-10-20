#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b3.list <<'L'
51 42 24 3 3 orbit 0 0 543 7
52 43 34 4 4 orbit 0 0 543 7
53 0 0 3 3 orbit 0 0 543 7
54 44 44 3 3 orbit 0 0 543 7
55 42 24 3 3 orbitT 0 0 543 7
56 0 0 3 3 orbitT 0 0 543 7
57 42 24 3 3 orbitN 0 0 543 7
58 0 0 0 0 orbitN 0 0 543 7
59 42 24 3 3 orbit 0 0 543 8
60 43 24 4 3 orbit 0 0 543 7
L
cat b3.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b3_s$0_$5_$9.txt 2>&1'
