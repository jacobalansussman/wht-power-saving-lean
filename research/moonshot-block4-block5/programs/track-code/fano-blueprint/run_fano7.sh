#!/bin/sh
# ms-fano: exact Fano-prefix searches at block 4 (one solved extra move per output). Block-3 sanity ran under bg.sh.
cd "$(dirname "$0")"
PY=python3
$PY -B fano7.py 4 PL 6 20 34 10 out_f7_m4_PL_k6.json       # control (8 moves)
$PY -B fano7.py 4 FN 5 400 35 10 out_f7_m4_FN_k5.json      # 2 in 7 target (mid gadget)
$PY -B fano7.py 4 FE 5 400 37 10 out_f7_m4_FE_k5.json      # 2 in 7 target (mirror end)
$PY -B fano7.py 4 FS 5 400 36 10 out_f7_m4_FS_k5.json      # 2 in 7 target (start gadget)
$PY -B fano7.py 4 FN 6 150 38 10 out_f7_m4_FN_k6.json      # 8-move count with the gadget (power check)
FANO_W=3 $PY -B fano7.py 4 FS 7 150 39 10 out_f7_m4_FS_W3_k7.json   # 3 in 10 target
