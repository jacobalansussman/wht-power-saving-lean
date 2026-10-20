#!/bin/sh
V=${PYTHON:-python3}   # a Python with python-sat (pip install python-sat)
cd "$(dirname "$0")"
for arg in "s2 9" "s2 10" "word 9" "word 10"; do
  echo "=== notW0 $arg $(date +%H:%M:%S)"
  perl -e 'alarm 1200; exec @ARGV' $V -I sat_notW0.py $arg 2>&1
  echo "exit $?"
done
