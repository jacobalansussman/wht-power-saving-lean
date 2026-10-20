#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ZTRIES=300 LMSTARTS=4 LMITERS=100
cat > b16.list <<'L'
709 42 24 3 3 orbitSNZ 4 3 3 8 300
801 42 24 3 3 orbitSNZ 4 3 3 8 300
802 42 24 3 3 orbitSNZ 4 3 3 8 300
803 42 24 3 3 orbitSNZ 4 3 3 8 300
804 42 24 3 3 orbitSNZ 4 3 3 8 300
805 42 24 3 3 orbitSNZ 4 3 3 8 300
L
cat b16.list | xargs -P 2 -L 1 bash -c 'timeout 600 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} > out_b16_s$0_$5_$9.txt 2>&1'
