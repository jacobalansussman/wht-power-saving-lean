#!/bin/bash
ulimit -v 12000000
python3 -B bf7c.py 4 61 20000 30000 > out_bf7c_s61.txt 2>&1 &
python3 -B bf7c.py 4 62 20000 30000 > out_bf7c_s62.txt 2>&1 &
wait
