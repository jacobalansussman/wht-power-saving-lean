#!/bin/sh
# F_2 decoders (XOR model) in the root-orthogonal flat family; each run capped at 25 min
V=${PYTHON:-python3}   # a Python with python-sat (pip install python-sat)
cd "$(dirname "$0")"
for prof in 6x8,7x2 6x9,8x1; do
  echo "=== xor s2 $prof $(date +%H:%M:%S)"
  perl -e 'alarm 1500; exec @ARGV' $V -I sat_rflats.py s2 $prof 2>&1
  echo "exit $?"
done
