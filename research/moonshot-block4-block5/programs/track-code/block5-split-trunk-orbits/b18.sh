#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ZTRIES=300 LMSTARTS=4 LMITERS=100
cat > b18.list <<'L'
1001 42 24 3 3 orbitSNZ4 0 0 543 7
1002 42 24 3 3 orbitSNZ4 0 0 543 7
1003 42 24 3 3 orbitSNZ4 0 0 543 7
1004 0 0 3 3 orbitSNZ4 0 0 543 7
1005 42 24 3 3 orbitSNZ4 0 0 543 8
1006 42 24 3 3 orbitSNZ4 0 0 543 8
1007 42 24 3 3 orbitSNZ4 0 0 543 8
1008 0 0 3 3 orbitSNZ4 0 0 543 8
1009 44 44 3 3 orbitSNZ4 0 0 543 7
1010 44 44 3 3 orbitSNZ4 0 0 543 8
L
cat b18.list | xargs -P 1 -L 1 bash -c 'timeout 400 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b18_s$0_$5_$9.txt 2>&1'
