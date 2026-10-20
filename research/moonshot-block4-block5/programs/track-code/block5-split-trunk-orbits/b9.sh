#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b9.list <<'L'
201 42 24 3 3 orbitCN 0 0 543 7
202 42 24 3 3 orbitCN 0 0 543 7
203 0 0 3 3 orbitCN 0 0 543 7
204 44 44 3 3 orbitCN 0 0 543 7
205 43 34 4 4 orbitCN 0 0 543 7
206 42 24 3 3 orbitCN 0 0 543 8
207 0 0 3 3 orbitCN 0 0 543 8
208 0 0 0 0 orbitCN 0 0 543 7
209 42 34 3 4 orbitCN 0 0 543 7
210 43 24 4 3 orbitCN 0 0 543 7
L
cat b9.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b9_s$0_$5_$9.txt 2>&1'
