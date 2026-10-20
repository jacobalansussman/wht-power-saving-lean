#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b8.list <<'L'
101 42 24 3 3 orbitC 0 0 543 7
102 42 24 3 3 orbitC 0 0 543 7
103 0 0 3 3 orbitC 0 0 543 7
104 44 44 3 3 orbitC 0 0 543 7
105 43 34 4 4 orbitC 0 0 543 7
106 42 24 3 3 orbitC 0 0 543 8
107 0 0 3 3 orbitC 0 0 543 8
108 0 0 0 0 orbitC 0 0 543 7
109 42 34 3 4 orbitC 0 0 543 7
110 43 24 4 3 orbitC 0 0 543 7
L
cat b8.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b8_s$0_$5_$9.txt 2>&1'
