#!/bin/sh
# x-pack: regenerate the Lean data of the best tier-1 variant (arithmetic + text generation only; NO Lean build here).
set -e
X=$(cd "$(dirname "$0")" && pwd)   # this folder (tools/community/repack)
: "${SHEA:?set SHEA to a clone of github.com/shea256/fourier-transform-below-nlogn at 5239484}"
: "${VERIFIER:?set VERIFIER to a clone of github.com/bigospeedrun/verifier}"
C=$SHEA/lean/complex588
REC=$VERIFIER/submissions/dft/nmamano-dft-record/Submission
cd $X
mkdir -p out/best
python3 -I tools/mk_price_json.py out/opt/c588x.schedule.json.gz out/best/c588x.json
for V in A-cubic B-quartic; do
  rm -rf lean-out/$V; mkdir -p lean-out/$V/Work/GCert/Data/Gen lean-out/$V/logs
  TK=tools; [ $V = B-quartic ] && TK=tools4
  python3 -I tools/generate_global_lean.py out/best/c588x.json --toolkit $TK --record $REC --out lean-out/$V | tee lean-out/$V/logs/generate.txt
  A=$(python3 -c "import json;print(json.load(open('lean-out/$V/logs/global-C588.json'))['kernel_numerator'])")
  D=$(python3 -c "import json;print(json.load(open('lean-out/$V/logs/global-C588.json'))['dft_numerator'])")
  python3 -I tools/patch_extra.py $C lean-out/$V 5009 $A $D
  if [ $V = B-quartic ]; then
    mkdir -p lean-out/$V/Work/BlockAccounting
    cp lean/Work/BlockAccounting/Moment4.lean lean-out/$V/Work/BlockAccounting/Moment4.lean
    python3 - lean-out/$V/Work/GCert/Data/RateC588x.lean <<'PY'
import sys
p=sys.argv[1]; t=open(p).read()
assert t.startswith("import Work.BlockAccounting.Moment\n")
open(p,'w').write(t.replace("import Work.BlockAccounting.Moment\n","import Work.BlockAccounting.Moment4\n",1))
PY
  fi
  # doc comment only: the generator's fixed text says 'cyclic schedule'; this table is an explicit edge-coloured one
  python3 - lean-out/$V/Work/GCert/Data/RateC588x.lean <<'PY'
import sys
p=sys.argv[1]; t=open(p).read()
a="and an explicit cyclic schedule. -/"
assert t.count(a)==1
open(p,'w').write(t.replace(a,"and an explicit 5,009-bank schedule (pad-optimal; x-pack). -/"))
PY
  (cd lean-out/$V && find Work -name '*.lean' | sort | xargs shasum -a 256 > logs/SHA256SUMS)
done
