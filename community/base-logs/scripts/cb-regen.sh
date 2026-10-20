#!/bin/bash
# c-base: the two Python steps of Shea's build_and_check.sh (provenance hashes; regenerate 63 Lean modules byte for byte).
P=/content/work/shea/lean/complex588
L=/content/work/logs/cb-build
mkdir -p $L
cd $P || exit 90
for s in check_provenance regenerate; do
  t=$(date +%s)
  echo "=== $s START $(date -u)"
  /usr/bin/time -v -o $L/py-$s.time python3 scripts/$s.py > $L/py-$s.log 2>&1
  rc=$?
  echo "=== $s END rc=$rc wall=$(( $(date +%s)-t ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/py-$s.time) last: $(tail -1 $L/py-$s.log)"
done
echo "git status lines after: $(git -C /content/work/shea status --short | wc -l)"
