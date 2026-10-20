#!/bin/bash
# t42 (resumed): pickle the d=4 and d=3 systems of union design s1401 (cap 4000: enough for balls of radius 2-3)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=3 LMITERS=80 CAP=4000 PROG=10000 MODPMAX=4e8 LMMAX=3000
DRV=jd_union ZTRIES=300 timeout 1500 python3 dumpsys.py 1401 42 24 3 3 orbitSNZ 0 0 543 8 > out_t42_s1401.txt 2>&1
