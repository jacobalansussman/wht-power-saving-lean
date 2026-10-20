#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ZTRIES=300 LMSTARTS=6 LMITERS=150
cat > b15a.list <<'L'
801 42 24 3 3 orbitSNZ 0 0 543 8
802 42 24 3 3 orbitSNZ 0 0 543 8
803 42 24 3 3 orbitSNZ 0 0 543 8
804 42 24 3 3 orbitSNZ 0 0 543 8
805 42 24 3 3 orbitSNZ 0 0 543 8
806 0 0 3 3 orbitSNZ 0 0 543 8
807 0 0 3 3 orbitSNZ 0 0 543 8
808 44 44 3 3 orbitSNZ 0 0 543 8
809 42 24 3 3 orbitSN 0 0 543 8
L
cat b15a.list | xargs -P 2 -L 1 bash -c 'timeout 600 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b15a_s$0_$5_$9.txt 2>&1'
