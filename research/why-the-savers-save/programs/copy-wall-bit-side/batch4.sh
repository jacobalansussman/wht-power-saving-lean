#!/bin/sh
V=${PYTHON:-python3}   # a Python with python-sat (pip install python-sat)
cd "$(dirname "$0")"
echo "=== sample s2 $(date +%H:%M:%S)"; perl -e 'alarm 300; exec @ARGV' $V -I sample_caseB.py s2 150 7 > out/sample_s2.out 2>&1; echo "exit $?"
echo "=== sample word $(date +%H:%M:%S)"; perl -e 'alarm 240; exec @ARGV' $V -I sample_caseB.py word 120 11 > out/sample_word.out 2>&1; echo "exit $?"
echo "=== xor rflats s2 6x8,7x2 $(date +%H:%M:%S)"; perl -e 'alarm 1500; exec @ARGV' $V -I sat_rflats.py s2 6x8,7x2 2>&1; echo "exit $?"
