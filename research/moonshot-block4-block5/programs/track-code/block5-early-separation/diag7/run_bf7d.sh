#!/bin/bash
ulimit -v 12000000
python3 -B bf7d.py 3 71 20000 20000 4 > out_bf7d_s71.txt 2>&1 &
python3 -B bf7d.py 3 72 20000 20000 4 > out_bf7d_s72.txt 2>&1 &
wait
