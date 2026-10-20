#!/bin/bash
# t35 (resumed): pickle the d=4 and d=3 systems of step-2 union designs s1401, s1402, same env as t12b
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=3 LMITERS=80 CAP=8000 PROG=10000 MODPMAX=4e8 LMMAX=3000
for s in 1401 1402; do
DRV=jd_union ZTRIES=300 timeout 900 python3 dumpsys.py $s 42 24 3 3 orbitSNZ 0 0 543 8 > out_t35_s$s.txt 2>&1
done
