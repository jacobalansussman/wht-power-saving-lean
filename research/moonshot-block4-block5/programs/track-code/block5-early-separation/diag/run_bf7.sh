#!/bin/bash
ulimit -v 12000000
python3 -B bf7.py 16 11 2 3 3000 > out_bf7_s11.txt 2>&1 &
python3 -B bf7.py 16 12 1 2 3000 > out_bf7_s12.txt 2>&1 &
python3 -B bf7.py 16 13 3 4 3000 > out_bf7_s13.txt 2>&1 &
wait
