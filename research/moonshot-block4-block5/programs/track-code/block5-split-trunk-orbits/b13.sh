#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b13.list <<'L'
601 42 24 3 3 orbitSN 0 0 543 7
602 42 24 3 3 orbitSN 0 0 543 7
603 42 24 3 3 orbitSN 0 0 543 7
604 0 0 3 3 orbitSN 0 0 543 7
605 44 44 3 3 orbitSN 0 0 543 7
606 42 24 3 3 orbitSN 0 0 543 8
607 42 24 3 3 orbitSN 0 0 543 8
608 42 24 3 3 orbitSNZ 0 0 543 7
609 42 24 3 3 orbitSNZ 0 0 543 8
610 0 0 3 3 orbitSNZ 0 0 543 8
611 0 0 0 0 orbitSN 0 0 543 7
612 42 24 3 3 orbitSNZ 0 0 543 7
L
cat b13.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b13_s$0_$5_$9.txt 2>&1'
