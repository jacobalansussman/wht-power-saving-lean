#!/bin/bash
export OMP_NUM_THREADS=1
time python3 rc_star5.py --noold > out_star5_noold.txt 2>&1
