#!/bin/bash
# Calibration of gen_mj / check_mj on small cases.  Usage: ./calibrate_depth.sh PYTHON_WITH_NUMPY
# Needs gen_m0, gen_mj, check_mj built next to it (see README) and dense_small.py.
cd "$(dirname "$0")"
PY=${1:-python3}
echo '1. With every bit a letter (lab = 2^m - 1) gen_mj is gen_m0: same rank and the same pivot rows and columns.'
for mk in '2 1' '3 4'; do
  set -- $mk; full=$(( (1 << $1) - 1 ))
  ./gen_m0 $1 $2 0 cal_a.txt >/dev/null 2>&1; ./gen_mj $1 $2 $full 0 cal_b.txt >/dev/null 2>&1
  if cmp -s <(tail -n +2 cal_a.txt) <(tail -n +2 cal_b.txt); then
    echo "   m=$1 k=$2: pivots identical ($(head -1 cal_b.txt | sed 's/.*size //') of them)"
  else echo "   m=$1 k=$2: DIFFERENT"; fi
done
echo "   check_mj on the block-3 one: $(./check_mj cal_b.txt 2>/dev/null | tail -1)"
echo
echo '2. The sum over all 2^m characters of rank M_alpha (gen_mj) against the rank of the whole block computed densely'
echo '   in floating point (dense_small.py, numpy SVD).'
for c in '3 4 1' '3 4 3' '3 4 7' '3 2 3' '4 2 1' '4 2 3' '4 2 7' '4 2 15'; do
  set -- $c; n=$((1 << $1)); tot=0; per=''
  for a in $(seq 0 $((n - 1))); do
    r=$(./gen_mj $1 $2 $3 $a /dev/null 2>/dev/null | sed 's/.*= //'); tot=$((tot + r)); per="$per $r"
  done
  echo "   gen_mj m=$1 k=$2 lab=$3: sum $tot (per character:$per)"
  echo "   $("$PY" -I dense_small.py $1 $2 $3)"
done
echo
echo '3. Block 4, k = 6, against hookG.py (data/block4-hook/out_rank2_blocks.txt: f >= 27,984 for b01 and 38,064 for b012).'
for lab in 3 7; do
  line=$(./gen_mj 4 6 $lab 0 /dev/null 2>/dev/null)
  echo "   $line  -> x 16 = $(( 16 * ${line##*=} ))"
done
