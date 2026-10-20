#!/bin/sh
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
M=${1:-35}
$VP -I sa2.py 11 $M 0.0 > sa2_11.log 2>&1 &
$VP -I sa2.py 12 $M 0.0 > sa2_12.log 2>&1 &
$VP -I sa2.py 13 $M 1.5 > sa2_13.log 2>&1 &
$VP -I sa2.py 14 $M 3.0 > sa2_14.log 2>&1 &
wait
