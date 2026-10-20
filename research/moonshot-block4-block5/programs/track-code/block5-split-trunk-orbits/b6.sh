#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b6.list <<'L'
71 42 24 3 3 orbitZ 0 0 543 7
72 0 0 3 3 orbitZ 0 0 543 7
73 44 44 3 3 orbitZ 0 0 543 7
74 42 24 3 3 orbitSZ 0 0 543 7
75 0 0 3 3 orbitSZ 0 0 543 7
76 44 44 3 3 orbitSZ 0 0 543 7
77 42 24 3 3 orbitZ 0 0 543 8
78 42 24 3 3 orbitSZ 0 0 543 8
79 42 24 3 3 orbitZ 0 0 543 7
80 42 24 3 3 orbitSZ 0 0 543 7
L
cat b6.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b6_s$0_$5_$9.txt 2>&1'
