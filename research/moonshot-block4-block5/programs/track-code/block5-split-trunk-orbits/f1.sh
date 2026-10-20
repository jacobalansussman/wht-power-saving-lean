#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for x in "61 42 24 3 3 orbitF 7" "64 44 44 3 3 orbitF 7" "68 44 44 3 3 orbitS 7" "104 44 44 3 3 orbitC 7" "65 42 24 3 3 orbitS 7" "77 42 24 3 3 orbitZ 8"; do
  set -- $x; timeout 600 python3 famcheck.py $x > out_f1_s$1.txt 2>&1
done
