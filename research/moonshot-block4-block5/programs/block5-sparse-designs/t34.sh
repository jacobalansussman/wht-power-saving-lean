#!/bin/bash
# t34 (resumed): pickle the d=4 and d=3 systems of step-2 design s1501 (zcover), same env as t24
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=3 LMITERS=80 CAP=6000 PROG=0 MODPMAX=4e8 LMMAX=3000
DRV=jd_zcover ZTRIES=300 ZC_MAX=4 timeout 900 python3 dumpsys.py 1501 42 24 3 3 orbitSNZ 0 0 543 8 > out_t34_s1501.txt 2>&1
