#!/bin/bash
export OMP_NUM_THREADS=1
( python3 anneal_gf2.py 3 2 1 20000 20 > out_n3_s2.txt; python3 anneal_gf2.py 3 1 2 20000 20 > out_n3_s1.txt ) &
( timeout 1500 python3 anneal_gf2.py 4 4 3 30000 10 > out_n4_s4.txt; timeout 2400 python3 anneal_gf2.py 4 3 4 30000 30 > out_n4_s3.txt ) &
wait
