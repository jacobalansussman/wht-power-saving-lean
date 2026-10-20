#!/bin/sh
# pf-logging (14 Oct 2026): FULL-ROW certificates (no row sampling) for the widened hook bounds:
# block 6, k = 3 (f = 4): seven weight classes alpha = 1^w; block 5, k = 4 (f = 5): alpha = 0.
# One process at a time (it ran on a Linux machine of Google's Colab service). Results: out_gen_full.txt, out_check_full.txt.
cd "$(dirname "$0")"
cc -O2 -o gen_big gen_big.c && cc -O2 -o check_big check_big.c || exit 1
date > out_run_full.txt
for w in 0 1 2 3 4 5 6; do
  a=$(( (1 << w) - 1 ))
  [ -f done_m6_a$a ] && continue
  echo "gen start $(date)" >> err_gen_full.txt; ./gen_big 6 3 $a 0 1 cert_m6_k3_a$a.full.txt >> out_gen_full.txt 2>> err_gen_full.txt
  ./check_big cert_m6_k3_a$a.full.txt >> out_check_full.txt 2>> err_check_full.txt
  date >> out_check_full.txt; touch done_m6_a$a
done
if [ ! -f done_m5_a0 ]; then
  echo "gen start $(date)" >> err_gen_full.txt; ./gen_big 5 4 0 0 1 cert_m5_k4_a0.full.txt >> out_gen_full.txt 2>> err_gen_full.txt
  ./check_big cert_m5_k4_a0.full.txt >> out_check_full.txt 2>> err_check_full.txt
  date >> out_check_full.txt; touch done_m5_a0
fi
date >> out_run_full.txt; echo ALLDONE >> out_run_full.txt
