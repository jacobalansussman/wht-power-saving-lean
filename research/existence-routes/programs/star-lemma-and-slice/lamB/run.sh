#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# (1) U_3 at the odd lambdas with other move orders (Singer order 1,2,4,3,6,7,5 and two others)
for ord in 1,2,4,3,6,7,5 7,6,5,4,3,2,1 3,5,6,7,1,2,4; do for lam in 1 7; do
  (python3 lamsolve.py 3 2 $ord $lam 25 2 > out_U3ord_${ord//,/}_l$lam.txt 2>&1 &); done; done
# (2) block 2, W = 2: every multiset of 4 and 5 moves from {1,2,3}, every lambda
for ms in 1,1,2,2 1,2,3,3 1,1,2,3 1,2,2,3 1,1,3,3 2,2,3,3 1,2,3,1,2 1,2,3,3,3 1,1,2,2,3 3,3,3,3; do for lam in 0 1 2 3; do
  (python3 lamsolve.py 2 2 $ms $lam 12 3 > out_b2_${ms//,/}_l$lam.txt 2>&1 &); done; done
sleep 5; while pgrep -f lamsolve.py > /dev/null; do sleep 5; done
cat out_*.txt
