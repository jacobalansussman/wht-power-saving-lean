#!/bin/sh
# run inside heavy.sh: 4 SA seeds in parallel (one slot), $1 minutes each
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
M=${1:-30}
$VP -I sa.py 1 $M 0.0 > sa_1.log 2>&1 &
$VP -I sa.py 2 $M 0.0 > sa_2.log 2>&1 &
$VP -I sa.py 3 $M 1.0 > sa_3.log 2>&1 &
$VP -I sa.py 4 $M 2.0 > sa_4.log 2>&1 &
wait
