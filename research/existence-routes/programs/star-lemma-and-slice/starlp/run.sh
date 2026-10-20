#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P=python3
# calibration first (feasible: plain 1 in 5; 1 in 4 from block 4 lifted? ; infeasible: uniform W=7, i.e. c = 1/7)
( $P rc_star5.py --fix=1:1,2:1,4:1,8:1,16:1 ; $P rc_star5.py --fix=uniform:1/7 ; $P rc_star5.py --fix=uniform:1/6 ; \
  $P rc_star5.py --fix=uniform:0.149 ; $P rc_star5.py --fix=uniform:0.151 ) > out_calib_star5.txt 2>&1 &
( $P rc_star5.py ; $P rc_star5.py --noold ) > out_star5.txt 2>&1 &
wait
