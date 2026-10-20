#!/bin/sh
# ms-fano: depth profile (delta_1, delta_3, delta_4) of the block-4 hooks V (x) Lambda^k V, k = 6, 7, 9, 10, 5
# (k = 8 done in run_delta3.sh).  Largest block 6435 x 6435 (~0.6 GB).  Question: does any k give
# delta_3 / delta_1 > 2.5 (then two depth-3 targets need 6 moves: no 2 in 7 ends privately for both arrays).
cd "$(dirname "$0")"
PY=python3
NS=../common/nullsp
for k in 7; do
  for spec in "d1 1" "d3 1,2,4" "d4 1,2,4,8"; do
    set -- $spec
    HOOKK=$k $PY -I hookk.py dump k${k}$1 $2 && $NS M_k${k}$1.bin null_k${k}$1.bin > /dev/null && HOOKK=$k $PY -I hookk.py load k${k}$1
  done
done
