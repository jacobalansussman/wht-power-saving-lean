#!/bin/bash
# pf-logging: commutative-slice certificates, small cases and 2-in-7 at m=4
cd "$(dirname "$0")"
VPY=${PYTHON:-python3}   # a Python 3 with numpy
for args in "4 1 4" "2 2 3" "3 2 5" "3 3 8" "4 2 7"; do
  set -- $args
  $VPY -I gen_matroid.py $1 $2 $3 cert_m$1_W$2_n$3.bin 2>&1 | grep -v -i warning | grep -v "viol ="
  /usr/bin/python3 -I check_matroid.py $1 $2 $3 cert_m$1_W$2_n$3.bin
done
