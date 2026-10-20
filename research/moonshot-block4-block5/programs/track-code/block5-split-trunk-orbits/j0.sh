#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 600 python3 jdesign.py 31 42 24 3 3 t2 2 2 543 > out_j0_t2.txt 2>&1 &
timeout 600 python3 jdesign.py 31 42 24 3 3 orbit 0 0 543 > out_j0_orbit.txt 2>&1 &
wait
