#!/bin/sh
# pf-logging: certificates for the widened-game hook bounds at block 6 (k=3, f=4: 63W/31) and block 5 (k=4, f=5: 62W/29)
cd "$(dirname "$0")"
date
for w in 0 1 2 3 4 5 6; do
  a=$(( (1 << w) - 1 ))
  ./gen_big 6 3 $a 8000 $((100 + w)) cert_m6_k3_a$a.txt >> out_gen_b56.txt 2>> err_gen_b56.txt
  ./check_big cert_m6_k3_a$a.txt >> out_check_b56.txt 2>> err_check_b56.txt
  date >> out_check_b56.txt
done
./gen_big 5 4 0 6500 7 cert_m5_k4_a0.txt >> out_gen_b56.txt 2>> err_gen_b56.txt
./check_big cert_m5_k4_a0.txt >> out_check_b56.txt 2>> err_check_b56.txt
date
