#!/bin/bash
ulimit -v 12000000
python3 -B bf8s.py 4 51 0 20000 8 > out_bf8s_s51.txt 2>&1 &
python3 -B bf8s.py 4 52 3 20000 8 > out_bf8s_s52.txt 2>&1 &
wait
