#!/bin/sh
# v-ms-fano: second-prime and uniformity checks of the block-4 hook depth values (~0.6 GB peak, a few min each)
cd "$(dirname "$0")"
PY=python3
NS=../../common/nullsp
run() { # prime k name dirs
  PRIME=$1 HOOKK=$2 $PY -I hk.py dump $3 $4 && $NS M_$3.bin null_$3.bin > /dev/null && PRIME=$1 HOOKK=$2 $PY -I hk.py load $3
}
run 33554341 8 p2k8d3 1,2,4
run 33554341 8 p2k8d1 1
run 33554393 8 p1k8d3b 3,5,9
run 33554341 8 p2k8d4 1,2,4,8
run 33554393 10 p1k10d3 1,2,4
run 33554393 9 p1k9d3 1,2,4
run 33554393 9 p1k9d1 1
