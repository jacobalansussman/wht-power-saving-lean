#!/bin/bash
# Re-create the A-fc fill: tolerant placement of the class solution on the words of #2979 (74 replicas).
# Usage: bash run2979.sh HUB3815 HUBMAIN OUT
#   HUB3815  a folder holding research/staged544-clock-ceiling of CrocSwap/integer-mult-bounds (tree e867b76f, = #3815 head 77ffbf28)
#   HUBMAIN  a checkout of CrocSwap/integer-mult-bounds at 12fa8c5c844be21bfefb38a19503fbd867a11244
#   OUT      a new output folder
# Needs Python 3 with numpy and highspy (the recorded run used Python 3.13, numpy 2.1.3, highspy 1.15.1).
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd)
D=$(cd "$H/../../../data/a-fc" && pwd)
HUB3815=$1; HUBMAIN=$2; OUT=$3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0 FC_COST=1,1 FC_WIDEN4=0
export FC_FALLBACK_LOG="$OUT.spare-roles.json"
python3 -B "$H/tools/build_fc2.py" 2979 "$HUB3815" "$HUBMAIN" "$HUBMAIN/research/rank7-joint-endframes/allocation/ORDER.json" "$OUT" \
  --C 24 --a 0.0017884 --from-sol "$D/class-solution/SOL-MIP0.json" --placer tolerant
