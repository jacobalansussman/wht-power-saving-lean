#!/bin/bash
# t11: every layer-3 "closest object" of w7-b5 / v-w7-b5 re-decided with the fixed peel2 + alpha3 (exact on V4)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=2 LMITERS=40 CAP=20000 PROG=10000 MODPMAX=4e8 LMMAX=1500
cat > t11.list <<'L'
1 1108 42 24 3 3 orbitCYZ 0 0 543 8
1 1107 42 24 3 3 orbitCYZ 0 0 543 8
1 1109 42 24 3 3 orbitCYZ 0 0 543 8
1 1101 42 24 3 3 orbitCY 0 0 543 8
1 1102 42 24 3 3 orbitCY 0 0 543 8
1 1103 42 24 3 3 orbitCY 0 0 543 8
1 1104 42 24 3 3 orbitCY 0 0 543 8
1 1301 42 24 3 3 orbitCYZ 0 0 543 8
1 1303 42 24 3 3 orbitCYZ 0 0 543 8
1 1306 42 24 3 3 orbitCYZ 0 0 543 8
1 1307 0 0 3 3 orbitCYZ 0 0 543 8
1 1308 0 0 3 3 orbitCY 0 0 543 8
1 706 42 24 3 3 orbitSN 0 0 543 8
1 707 42 24 3 3 orbitSN 0 0 543 8
1 709 42 24 3 3 orbitSNZ 0 0 543 8
1 710 0 0 3 3 orbitSNZ 0 0 543 8
300 808 44 44 3 3 orbitSNZ 0 0 543 8
L
while read zt s rest; do
ZTRIES=$zt timeout 600 python3 alpha3.py $s $rest > out_t11_s$s.txt 2>&1
done < t11.list
