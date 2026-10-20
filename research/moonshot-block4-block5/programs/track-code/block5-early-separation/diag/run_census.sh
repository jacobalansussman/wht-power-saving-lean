#!/bin/bash
ulimit -v 12000000
for m in 3 4 5; do python3 -B census.py $m 12 3000 $m > out_census_m$m.txt 2>&1 & done
wait
