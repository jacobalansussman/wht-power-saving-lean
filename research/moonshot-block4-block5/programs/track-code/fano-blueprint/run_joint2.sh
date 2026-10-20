#!/bin/sh
# ms-fano: exact end-gadget completion (two SHARED solved moves) after sampled Fano prefixes.
cd "$(dirname "$0")"
PY=python3
$PY -B joint2.py 3 PL 4 10 21 3 out_j2_m3_PL_k4.json     # control: 6-move words
$PY -B joint2.py 3 FN 3 40 22 3 out_j2_m3_FN_k3.json     # block 3, 5 moves: must be none
$PY -B joint2.py 4 PL 6 10 23 3 out_j2_m4_PL_k6.json     # control: 8-move words
$PY -B joint2.py 4 FN 5 150 24 3 out_j2_m4_FN_k5.json    # 2 in 7 (mid gadget + solved shared end)
$PY -B joint2.py 4 FS 5 120 25 3 out_j2_m4_FS_k5.json    # 2 in 7 (start gadget + solved shared end)
$PY -B joint2.py 4 FN 6 60 26 3 out_j2_m4_FN_k6.json     # 8 moves with the gadget (power check)
