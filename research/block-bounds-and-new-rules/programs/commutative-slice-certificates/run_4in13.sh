#!/bin/bash
# pf-logging, 14 Oct 2026: commutative slice at block 4, 4 arrays in 13 moves (20,058,300 multisets).
cd "$(dirname "$0")"
VPY=${PYTHON:-python3}   # a Python 3 with numpy
date
/usr/bin/time -l $VPY -I gen_matroid.py 4 4 13 cert_m4_W4_n13.bin 2>&1 | grep -v -i warning | grep -v "viol =" | grep -v "^ \+[0-9]* \+[a-z]"
date
/usr/bin/time -l /usr/bin/python3 -I check_matroid.py 4 4 13 cert_m4_W4_n13.bin 2>&1 | grep -E "m=|exclusion|RESULT|maximum resident|real"
date
