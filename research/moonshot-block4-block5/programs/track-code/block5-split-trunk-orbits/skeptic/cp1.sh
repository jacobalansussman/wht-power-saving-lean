#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd cyc
cat > cp1.list <<'L'
1201 42 24 3 3 t2 0 0 54 8 0 s2b,moveA2
1202 42 24 3 3 t2 0 0 54 8 0 s2b,moveA2
1203 42 24 3 3 t2 0 0 54 7 0 s2b,moveA2
1204 42 24 3 3 t2 0 0 54 7 0 s2b,moveA2
1205 42 24 3 3 t2 0 0 54 8 0 s2b,fixNb
1206 42 24 3 3 t2 0 0 54 8 0 s2b
L
cat cp1.list | xargs -P 2 -L 1 bash -c 'timeout 400 python3 jdesign.py $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} > out_cp1_s$0_$5_$9.txt 2>&1'
