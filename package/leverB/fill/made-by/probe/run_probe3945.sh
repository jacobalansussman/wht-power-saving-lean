#!/bin/sh
# o-model probe P3 (C1): #3815's Hall LP with a head bonus on the child-8/10 empty pattern (scale $1), #3815's build_fill.py
# (patched only to import the patched hall_lp; same directory), owner_audit.py of #3734 on the result, then #3945's offset
# step (offsets3945.py) with #3945's final widths, then head count. NAME=$2.
set -e
B=/content/work/reopt-o-model; PY=$B/v3815/bin/python; SCALE=$1; NAME=$2
R=$B/c3815h-$NAME/research; rm -rf $B/c3815h-$NAME; mkdir -p $R/clock-heads/construction
cp $B/hub-3815/research/staged544-clock-ceiling/construction/build_fill.py $R/clock-heads/construction/
$PY -I -B $B/probe/mkpatch3815.py $B/hub-3815/research/staged544-clock-ceiling/construction/hall_lp.py $R/clock-heads/construction/hall_lp.py $SCALE
ln -s $B/hub-3815/research/staged544-producer-return $R/staged544-producer-return
cd $B/c3815h-$NAME
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
/usr/bin/time -v nice -n 19 $PY -B $R/clock-heads/construction/build_fill.py $B/hub-3815/research/staged544-clock-ceiling/inputs/ORDER.json --out $B/c3815h-$NAME/out --distinct 40
nice -n 19 $PY -I -B $B/tools/cmpfill3815.py $B/hub-3815/research/staged544-producer-return/inputs/NATIVE-base.json.gz $B/hub-3815/research/staged544-producer-return/inputs/TEMPLATE.json.gz $B/c3815h-$NAME/out/inputs/FILL.json.gz $B/hub-3815/research/staged544-clock-ceiling/inputs/FILL.json.gz
nice -n 19 $PY -I -B $B/probe/offsets3945.py $B/hub-3945/research/staged544-terminal-gap-fill/owner $B/c3815h-$NAME/out/inputs/FILL.json.gz $B/c3815h-$NAME/final3945
echo DONE_P3 $NAME
