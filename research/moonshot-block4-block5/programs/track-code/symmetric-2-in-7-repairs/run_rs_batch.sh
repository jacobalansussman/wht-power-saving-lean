#!/bin/bash
# run_rs_batch.sh H W shape seeds... : relaxation-guided search (rs2.py, 300 steps) for several seeds, start 0 and 1,
# then projective normalisation + full test of every final skeleton whose relaxation residual < 1e-6
cd "$(dirname "$0")"
H=$1; W=$2; SH=$3; shift 3
for sd in "$@"; do for st in 0 1; do
  P=rsb_${H}_${sd}_${st}
  timeout 900 python3 rs2.py 4 $W $H $SH 2 $sd $st 300 $P > $P.log 2>&1
  python3 - <<PY >> $P.log 2>&1
import pickle
try:
    d = pickle.load(open("${P}_final.pkl", "rb")); print("FINALRES", d["res"])
except Exception as e:
    print("nofinal", e)
PY
  R=$(grep FINALRES $P.log | awk '{print $2}')
  if python3 -c "import sys; sys.exit(0 if float('${R:-1}') < 1e-6 else 1)"; then
    timeout 900 python3 normalize.py ${P}_final.pkl 1e-3 4 6000 >> $P.log 2>&1
  fi
done; done
echo BATCHDONE
