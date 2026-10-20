#!/bin/sh
# run_fg.sh "KIND W k n seed" ... : fano_gate runs, one file each
cd "$(dirname "$0")"
for a in "$@"; do
  set -- $a
  python3 -B fano_gate.py $1 $2 $3 $4 $5 > out_fg_$1_$2_$3_s$5.txt 2>&1
  tail -1 out_fg_$1_$2_$3_s$5.txt
done
