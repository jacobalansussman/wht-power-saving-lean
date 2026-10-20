#!/bin/sh
# o-model probe (C1): #2979's T37 flow model with head variables. P1: all child-8/10 bins empty (feasibility, no hint).
# P2: maximize heads, hinted with the published T37 flow. Then double.py -> T74 fill, owner engine check, head count.
set -e
B=/content/work/reopt-o-model; H=$B/hub-main/research; P=$H/rank7-joint-endframes; C=$B/c2979; PY=$B/v2979/bin/python
cd $C
$PY -I -B $B/probe/mkpatch.py code/fill17.py code/fill17_heads.py
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
MODE=$1; NAME=$2; shift 2
rm -rf runs/$NAME-t37 runs/$NAME-t74
/usr/bin/time -v nice -n 19 $PY -B code/fill17_heads.py --work runs/$NAME-t37 --chronology pr2833-t37 --providers $P/allocation/PROVIDERS.json \
  --mix-file $P/construction/allocation/T37-FACE-MIX.json --buffers 91 --seed 1 --workers 1 --catalogue tight $MODE "$@" \
  --method 'o-model probe: #2979 fill17.py + head variables (fill17_heads.py) | doubled by double.py'
/usr/bin/time -v nice -n 19 $PY -B code/double.py runs/$NAME-t37 runs/$NAME-t74 --providers $P/allocation/PROVIDERS.json --buffers 91
$PY -I -B $B/tools/groupstats.py runs/$NAME-t74/FILL.json.gz $P/allocation/ORDER.json $P/allocation/DISPATCH.json $P/allocation/PROVIDERS.json $P/allocation/vendor/TEMPLATE.json.gz > runs/$NAME-t74/GROUPSTATS.json
grep -A3 -E '"8/10|"12/13|total_heads' runs/$NAME-t74/GROUPSTATS.json
sha256sum runs/$NAME-t74/FILL.json.gz
echo DONE_PROBE $NAME
