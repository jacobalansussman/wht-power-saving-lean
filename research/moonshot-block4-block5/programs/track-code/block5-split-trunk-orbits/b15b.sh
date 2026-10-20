#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ZTRIES=300 LMSTARTS=6 LMITERS=150
cat > b15b.list <<'L'
901 42 24 3 3 orbitSNZ 0 0 543 7
902 42 24 3 3 orbitSNZ 0 0 543 7
903 42 24 3 3 orbitSNZ 0 0 543 7
904 42 24 3 3 orbitSNZ 0 0 543 7
905 42 24 3 3 orbitSNZ 0 0 543 7
906 0 0 3 3 orbitSNZ 0 0 543 7
907 0 0 3 3 orbitSNZ 0 0 543 7
908 44 44 3 3 orbitSNZ 0 0 543 7
909 42 24 3 3 orbitSN 0 0 543 7
L
cat b15b.list | xargs -P 2 -L 1 bash -c 'timeout 600 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b15b_s$0_$5_$9.txt 2>&1'
