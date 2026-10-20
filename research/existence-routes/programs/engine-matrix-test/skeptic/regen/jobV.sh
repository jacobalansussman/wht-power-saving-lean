#!/bin/bash
# skeptic v-e3-engineE: regenerate e3-engineE's block-4 realised fake (their hecke_vp.py, same args and seed as their
# jobF3.sh line "4 2 5 8 3", plus a save file) and check it with my own frame module (v_vpcheck.py).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NOPERM=1
{ timeout 1500 python3 -W ignore hecke_vp.py 4 2 5 8 3 0 vp_4_2_5.npz; timeout 900 python3 -W ignore v_vpcheck.py vp_4_2_5.npz; echo DONE; } > out_jobV.txt 2>&1
