#!/bin/sh
# matching (Hall) relaxation, field-free; each run capped at 20 min
V=${PYTHON:-python3}   # a Python with python-sat (pip install python-sat)
cd "$(dirname "$0")"
for prof in 6x8,7x2 6x9,8x1 7x8; do
  echo "=== match s2 $prof $(date +%H:%M:%S)"
  perl -e 'alarm 1200; exec @ARGV' $V -I match_sat.py s2 $prof 2>&1
  echo "exit $?"
done
