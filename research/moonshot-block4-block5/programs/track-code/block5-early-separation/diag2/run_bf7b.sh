#!/bin/bash
ulimit -v 12000000
python3 -B bf7b.py 6 21 2 3 20000 > out_bf7b_s21.txt 2>&1 &
python3 -B bf7b.py 6 22 1 2 20000 > out_bf7b_s22.txt 2>&1 &
python3 -B bf7b.py 6 23 3 4 20000 > out_bf7b_s23.txt 2>&1 &
wait
