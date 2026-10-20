#!/bin/bash
# t36 (resumed): calibration input for v4gb.py: pickle s801's d=4 system (binomial; alpha3 knows its V4: 8 families)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONHASHSEED=0 LMSTARTS=1 LMITERS=60 CAP=20000 ZTRIES=300 MODPMAX=4e8 LMMAX=3000
DRV=alpha3 timeout 600 python3 dumpsys.py 801 42 24 3 3 orbitSNZ 0 0 54 8 > out_t36_s801.txt 2>&1
