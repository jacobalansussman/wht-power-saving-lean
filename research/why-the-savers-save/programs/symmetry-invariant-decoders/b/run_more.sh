#!/bin/sh
# t-km: further cyclic-group KM runs (flats_cache.npz).  Each line under heavy.sh.
H=checks/wht42/common/heavy.sh
VP=../venv/bin/python
$H $VP -I kmg.py 8 8 6 63 > kmg_8_6.out 2>&1
SMALLTYPES=1x4,1x7,1x8,7x1 MAXPIECES=4 $H $VP -I kmg.py 7 7 5 63 > kmg_7_5.out 2>&1
SECOND=4x7 $H $VP -I kmg.py 4 4 7 63 > kmg_4_7_4x7.out 2>&1
SMALLTYPES=1x6,1x7,2x3,2x4,2x5,2x6,2x7,4x2,4x3,8x1 MAXPIECES=3 $H $VP -I kmg.py 8 8 5 63 > kmg_8_5.out 2>&1
SECOND=4x7 MAXPIECES=3 $H $VP -I kmg.py 4 4 6 63 > kmg_4_6_4x7.out 2>&1
SMALLTYPES=1x4,1x7,1x8,7x1 MAXPIECES=5 $H $VP -I kmg.py 7 7 4 63 > kmg_7_4.out 2>&1
SECOND=3x7 SMALLTYPES=3x7 MAXPIECES=1 $H $VP -I kmg.py 3 3 7 63 > kmg_3_7_3x7_3x7.out 2>&1
