#!/bin/sh
V=${PYTHON:-python3}   # a Python with python-sat (pip install python-sat)
cd "$(dirname "$0")"
for arg in "s2 6x2 --caseB" "s2 6x8 --caseB" "s2 6x8"; do
  echo "=== sat_light $arg $(date +%H:%M:%S)"; perl -e 'alarm 1200; exec @ARGV' $V -I sat_light.py $arg 2>&1; echo "exit $?"
done
