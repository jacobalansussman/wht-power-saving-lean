#!/bin/bash
# k-lower: build the new modules ONE AT A TIME in MODULES.txt order (LEAN_NUM_THREADS=2), waiting while 2 or more Lean
# processes run on the machine; one log per module under /content/work/logs/k-lower/; then the checks.
set -uo pipefail
. /content/work/env.sh
PKG=/content/work/k-lower-pkg-u8x
T=/content/work/k-lower/lean/complex588
L=/content/work/logs/k-lower-u8x
mkdir -p $L
cd $T
echo "BUILD START $(date -u)"
n=0
while read m; do
  [ -z "$m" ] && continue
  n=$((n+1))
  while [ "$(pgrep -c -x lean)" -ge 2 ]; do sleep 5; done
  log=$L/$(printf %02d $n)-$m.log
  LEAN_NUM_THREADS=2 /usr/bin/time -v lake build $m > $log 2>&1
  rc=$?
  t=$(grep -E "Elapsed" $log | tail -1 | awk '{print $NF}'); r=$(grep -E "Maximum resident" $log | awk '{print $NF}')
  echo "$(printf %02d $n) rc=$rc $m wall=$t maxrss_kb=$r"
  if [ $rc -ne 0 ]; then echo "FAILED at $m (log $log)"; grep -E "error" $log | head -20; exit 1; fi
done < $PKG/MODULES.txt
echo "== axioms"
LEAN_NUM_THREADS=2 lake env lean Work/FourierU8X/Axioms.lean 2>&1 | tee $L/axioms.log
echo "== zero / exit / glue axioms"
cat > /tmp/k-lower-ax.lean <<'LEAN'
import Work.GCert.Data.U8XGlue401
#print axioms GXD.U8X.exitRC
#print axioms GXD.U8X.ginv_cost
#print axioms GXD.U8X.fitsS
#print axioms GXD.U8X.price
#print axioms GXD.U8X.validX
LEAN
cp /tmp/k-lower-ax.lean Work/KLow/AxCheck.lean
LEAN_NUM_THREADS=2 lake env lean Work/KLow/AxCheck.lean 2>&1 | tee $L/axioms-glue.log
rm -f Work/KLow/AxCheck.lean
echo "== sorry grep"
grep -rnw sorry --include=*.lean Work/KLow Work/GCert/Data/*U8X* Work/GCert/Data/Gen/U8X* Work/FourierU8X | grep -v '`sorry`' | tee $L/sorry-grep.log
echo "== statement check (Compare.lean), positive then the one-unit-stronger negative"
while [ "$(pgrep -c -x lean)" -ge 1 ]; do sleep 5; done
LEAN_NUM_THREADS=2 lake build Work.FourierU8X.UniformFourierChallenge Work.FourierU8XNeg.UniformFourierChallenge > $L/challenges.log 2>&1; echo "challenges rc=$?"
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.FourierU8X.UniformFourierChallenge Work.FourierU8X.Main \
  OAI.PowerSaving.transform_mainU8X OAI.PowerSaving.convolution_mainU8X > $L/compare.log 2>&1; echo "compare rc=$?"; grep RESULT $L/compare.log
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.FourierU8XNeg.UniformFourierChallenge Work.FourierU8X.Main \
  OAI.PowerSaving.transform_mainU8X OAI.PowerSaving.convolution_mainU8X > $L/compare-neg.log 2>&1; echo "compare-neg rc=$? (must FAIL)"; grep RESULT $L/compare-neg.log
echo "BUILD END $(date -u)"
