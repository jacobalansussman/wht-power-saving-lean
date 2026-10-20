#!/bin/bash
ulimit -v 12000000
python3 -B bf7b.py 6 31 2 3 20000 > out_bf7b_s31.txt 2>&1 &
python3 -B bf7b.py 6 32 1 2 20000 > out_bf7b_s32.txt 2>&1 &
wait
