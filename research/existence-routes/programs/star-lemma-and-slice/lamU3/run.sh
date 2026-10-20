#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# U_3 = 2 in 7 at every lambda; plain 2 in 6 (control, exists); 2 in 5 = {1,2,4,3,5}? (control, impossible)
for lam in 0 1 2 3 4 5 6 7; do (python3 lamsolve.py 3 2 1,2,3,4,5,6,7 $lam 30 1 > out_U3_l$lam.txt 2>&1 &) ; done
for lam in 0 1 3 7; do (python3 lamsolve.py 3 2 1,1,2,2,4,4 $lam 10 1 > out_plain26_l$lam.txt 2>&1 &); done
sleep 2; while pgrep -f lamsolve.py > /dev/null; do sleep 5; done
cat out_*.txt
