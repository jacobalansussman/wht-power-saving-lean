#!/bin/sh
# ms-fano: delta_1, delta_2 (calibration) and delta_3 of the f = 9 block-4 hook. ~0.6 GB, a few minutes each.
cd "$(dirname "$0")"
PY=python3
NS=../common/nullsp
for spec in "d3 1,2,4" "d1 1" "d2 1,2"; do
  set -- $spec
  $PY -I delta3.py dump $1 $2 && $NS M_$1.bin null_$1.bin > /dev/null && $PY -I delta3.py load $1
done
