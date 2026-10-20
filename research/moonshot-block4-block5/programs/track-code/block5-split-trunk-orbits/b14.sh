#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cat > b14.list <<'L'
701 42 24 3 3 orbitSN 0 0 543 7
702 42 24 3 3 orbitSN 0 0 543 7
703 42 24 3 3 orbitSN 0 0 543 7
704 0 0 3 3 orbitSN 0 0 543 7
705 44 44 3 3 orbitSN 0 0 543 7
706 42 24 3 3 orbitSN 0 0 543 8
707 42 24 3 3 orbitSN 0 0 543 8
708 42 24 3 3 orbitSNZ 0 0 543 7
709 42 24 3 3 orbitSNZ 0 0 543 8
710 0 0 3 3 orbitSNZ 0 0 543 8
711 0 0 0 0 orbitSN 0 0 543 7
712 42 24 3 3 orbitSNZ 0 0 543 7
L
cat b14.list | xargs -P 2 -L 1 bash -c 'timeout 900 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 > out_b14_s$0_$5_$9.txt 2>&1'
