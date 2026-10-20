#!/bin/bash
# t46 (resumed): msolve on LEAF-PEELED ball tests (exact reduction).  s801 calibration (k2 must be [-1]; b4 positive-dim),
# then s1501 k2, k1, b4.  [-1] = no solution over the algebraic closure of F_p = EXCLUDED.
for f in L801_k2 L801_b4 L1501_k2 L1501_b4 L1501_k1; do
  timeout 1200 msolve -t 1 -v 1 -f $f.ms -o out_t46_$f.res > out_t46_$f.txt 2>&1
  echo "$f: $(head -c 60 out_t46_$f.res | tr "\n" " ") (${SECONDS}s cumulative)" >> out_t46_summary.txt
done
