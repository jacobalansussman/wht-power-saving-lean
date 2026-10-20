#!/bin/sh
# k-word: reproduce the hub's U8 canonical word (reproduce_barch.py, RECIPE-U8) in /content/work/k-word
K=/content/work/k-word; HUB=$K/hub; L=/content/work/logs/k-word
SRC=$HUB/research/complex-k492-owner-repair/vendor/multicold
cd $HUB
export PYTHONDONTWRITEBYTECODE=1 PYTHONUTF8=1 PYTHONINTMAXSTRDIGITS=0 MULTICOLD_SOURCE=$SRC GAIN_PR564_SOURCE=$SRC
unset PYTHONOPTIMIZE
echo "START $(date -u +%FT%TZ) hub $(git rev-parse HEAD)" > $L/u8.log
/usr/bin/time -v python3 -B $SRC/code/run_portable.py $HUB/research/complex-t21-owners/vendor/pr870/construction/barch/reproduce_barch.py \
  --main $HUB --pr635 $HUB/research/complex-sw4-frame8/vendor/pr635 --jx-dir $K/jx799 --recipe U8 --out $K/u8-native >> $L/u8.log 2>&1
echo "EXIT $? $(date -u +%FT%TZ)" >> $L/u8.log
