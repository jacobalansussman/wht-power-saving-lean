#!/bin/bash
# fc-reopt round 2: tolerant placement of a given class solution on #2979's line (T74). args: NAME SOL A [extra args]
B=/content/work/fc-reopt; NAME=$1; SOL=$2; A=$3; shift 3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 FC_COST=1,1 FC_WIDEN4=0
mkdir -p $B/r2/runs; export FC_FALLBACK_LOG=$B/r2/runs/$NAME.fallback.json
ORDER=$B/hub-main/research/rank7-joint-endframes/allocation/ORDER.json
cd $B && /usr/bin/time -v nice -n 10 python3 -B $B/r2/tools/build_fc2.py 2979 $B/hub-3815 $B/hub-main $ORDER $B/r2/runs/$NAME --C 24 --a $A --from-sol $SOL --placer tolerant "$@"
echo EXIT $?
