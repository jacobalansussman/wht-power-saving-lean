#!/bin/sh
# Reproduction of the bank re-packing of community/c588-repack (arithmetic + text generation; python3 -I; inputs read-only).  Lean: see the end.
set -e
X=$(cd "$(dirname "$0")" && pwd)   # this folder (tools/community/repack)
: "${SHEA:?set SHEA to a clone of github.com/shea256/fourier-transform-below-nlogn at 5239484}"
: "${VERIFIER:?set VERIFIER to a clone of github.com/bigospeedrun/verifier}"
C=$SHEA/lean/complex588; REC=$VERIFIER/submissions/dft/nmamano-dft-record/Submission
cd $X
# 0. CONTROL: Shea's own table read back from his Lean -> W=12208, 5,040 banks, kernel 12331686, DFT 12331685;
#    and his generator reproduces all 35 of his generated Lean files byte for byte from it
python3 -I tools/shea_sched_from_lean.py $C/Work/GCert/Data/Gen out/control/shea.schedule.json.gz
python3 -I tools/price_sched.py out/control/price-shea.json out/control/shea.schedule.json.gz
python3 -I tools/mk_price_json.py out/control/shea.schedule.json.gz out/control/shea.json
O=$(mktemp -d); mkdir -p $O/Work/GCert/Data/Gen $O/logs
python3 -I tools/generate_global_lean.py out/control/shea.json --toolkit tools --record $REC --out $O
(cd $O && for f in $(find Work -name '*.lean'); do cmp -s $f $C/$f || echo "DIFF $f"; done; echo "generator control: $(find Work -name '*.lean' | wc -l) files compared")
# 1. optimal packing: the integer programme (HiGHS in ./venv: python3 -m venv venv; venv/bin/pip install highspy) and the table
#    (10 min time limit; it finds the optimum, the proof of optimality is step 2)
venv/bin/python -I tools/pack_opt.py $C/Work/GCert/Data/Gen/C588Sched.lean out/opt/c588x.schedule.json.gz
python3 -I tools/price_sched.py out/opt/price-opt.json out/opt/c588x.schedule.json.gz
python3 -I tools/schedk_emulate.py out/opt/c588x.schedule.json.gz 5009
# 2. optimality proof of the pad invoice (invariant + exact enumeration), and the 4,999-bank 'ceiling' is infeasible
python3 -I tools/pad_optimal.py
venv/bin/python -I tools/pack_ilp.py $C/Work/GCert/Data/Gen/C588Sched.lean /dev/null || true   # 'kInfeasible' = no 4,999-bank table
# 3. rate certificate: third-order (term_le) vs fourth-order (term_le4) grid limits
python3 -I tools/rate_gap.py out/control/shea.schedule.json.gz out/opt/c588x.schedule.json.gz
# 4. Lean data for the Colab: lean-out/A-cubic (data swap only) and lean-out/B-quartic (+ Moment4.lean)
./MAKE-LEAN-OUT.sh
# Lean done on this Mac (light, via locked.sh): lake build Work.BlockAccounting.Moment4 ; lake build Work.XPackTest.Term4Test
