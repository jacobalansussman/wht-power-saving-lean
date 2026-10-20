#!/bin/bash
cd "$(dirname "$0")"
VPY=${PYTHON:-python3}   # a Python 3 with numpy
/usr/bin/time -l $VPY -I gen_matroid.py 4 3 10 cert_m4_W3_n10.bin 2>&1 | grep -v -i warning | grep -v "viol =" | grep -v "^ \+[0-9]* \+[a-z]" 
/usr/bin/time -l /usr/bin/python3 -I check_matroid.py 4 3 10 cert_m4_W3_n10.bin 2>&1 | grep -E "m=|exclusion|RESULT|maximum resident|real"
