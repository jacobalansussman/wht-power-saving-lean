#!/bin/bash
# bB-data: p-wgen's two #3945 provider words (ENG, fill-independent) and x-inst's generic dispatch Disp2979W, built in
# my tree NEXT TO phase A of bb-build.sh (disjoint modules; the shared dependencies are all built by the pre-build).
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
cd $T; L=$LOGROOT/words; mkdir -p $L
wl 14000
step words env LEAN_NUM_THREADS=${WTHREADS:-2} lake build Work.PWGen.S5Pair Work.XHostT.HFW2979.Disp2979W
rc=$?; grep -E "✖|error:" $L/words.log | head; tail -2 $L/words.log; echo "WORDS rc=$rc"
