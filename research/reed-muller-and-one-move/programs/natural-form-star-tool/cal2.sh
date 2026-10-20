#!/bin/bash
# cal2.sh: run starrec2 on "m ks [label]" lines from stdin (test forms = every class with r = 0 for these targets)
while read -r m ks label; do
  [ -z "$m" ] && continue
  vm=$(${PYTHON:-python3} vmask.py tables$m.bin $ks)
  r=$(echo "$ks" | VMASK=$vm ROUNDS=${ROUNDS:-3} ./starrec2 tables$m.bin list 2>/dev/null | grep "^S" | cut -c1-200)
  echo "m=$m [$label] $r"
done
