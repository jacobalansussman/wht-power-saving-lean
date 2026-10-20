#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd cyc
cat > cy1.list <<'L'
1101 42 24 3 3 orbitCY 0 0 543 8
1102 42 24 3 3 orbitCY 0 0 543 8
1103 42 24 3 3 orbitCY 0 0 543 8
1104 42 24 3 3 orbitCY 0 0 543 8
1105 42 24 3 3 orbitCY 0 0 543 7
1106 42 24 3 3 orbitCY 0 0 543 7
1107 42 24 3 3 orbitCYZ 0 0 543 8
1108 42 24 3 3 orbitCYZ 0 0 543 8
1109 42 24 3 3 orbitCYZ 0 0 543 8
1110 42 24 3 3 orbitCYZ 0 0 543 7
L
cat cy1.list | xargs -P 2 -L 1 bash -c 'timeout 600 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_cy1_s$0_$5_$9.txt 2>&1'
