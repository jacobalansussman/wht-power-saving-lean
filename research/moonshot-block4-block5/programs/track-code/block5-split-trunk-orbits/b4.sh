#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b4.list <<'L'
61 42 24 3 3 orbit 0 0 543 7
62 43 34 4 4 orbit 0 0 543 7
63 0 0 3 3 orbit 0 0 543 7
64 44 44 3 3 orbit 0 0 543 7
65 42 24 3 3 orbitS 0 0 543 7
66 43 34 4 4 orbitS 0 0 543 7
67 0 0 3 3 orbitS 0 0 543 7
68 44 44 3 3 orbitS 0 0 543 7
69 42 24 3 3 orbit 0 0 543 8
70 42 24 3 3 orbitS 0 0 543 8
L
cat b4.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b4_s$0_$5_$9.txt 2>&1'
