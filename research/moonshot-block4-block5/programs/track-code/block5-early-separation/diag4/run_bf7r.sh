#!/bin/bash
ulimit -v 12000000
python3 -B bf7r.py 4 41 2 20000 3 > out_bf7r_s41.txt 2>&1 &
python3 -B bf7r.py 4 42 1 20000 3 > out_bf7r_s42.txt 2>&1 &
wait
