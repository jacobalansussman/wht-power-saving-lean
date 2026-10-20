#!/bin/bash
export OMP_NUM_THREADS=1
python3 vorb.py 11 4 20000 > out_vorb_11.txt 2>&1 &
python3 vorb.py 12 4 20000 > out_vorb_12.txt 2>&1 &
wait
