#!/bin/sh
cd "$(dirname "$0")"
VP=${PYTHON:-python3}
M=${1:-40}
$VP -I sa3.py 21 $M 0.0 > sa3_21.log 2>&1 &
$VP -I sa3.py 22 $M 2.0 > sa3_22.log 2>&1 &
wait
