#!/bin/bash
# k-combine: Lake build of U8XM in my own tree (traces written, so the comparator need not rebuild), one Lake job at a
# time (LEAN_NUM_THREADS=1), started only when fewer than 2 heavy Lean processes run; then the checks; then the
# official comparator (positive, then negative).
. /content/work/env.sh
P=/content/work/kcomb-pkg/k-combine-u8xm
T=/content/work/k-combine/lean/complex588
L=/content/work/logs/k-combine
N=U8XM
cd $T
TG="Work.FourierU8XM.Main Work.FourierU8XM.UniformFourierChallenge Work.FourierU8XMNeg.UniformFourierChallenge Work.GCert.Data.U8XMControls Work.FourierU8XM.Axioms"
w1() { while [ "$(ps -C lean -o rss= | awk "\$1>2000000" | wc -l)" -ge 2 ]; do sleep 5; done; }
w0() { while [ "$(ps -C lean -o rss= | awk "\$1>2000000" | wc -l)" -ge 1 ]; do sleep 5; done; }
w1
echo "LAKE BUILD START $(date -u)" > $L/lake.out
LEAN_NUM_THREADS=1 /usr/bin/time -v lake build $TG >> $L/lake.out 2>&1
rc=$?
echo "LAKE BUILD END rc=$rc $(date -u)" >> $L/lake.out
[ $rc -ne 0 ] && { echo "RUN2 STOP (build failed)" >> $L/lake.out; exit 1; }
{
echo "== axioms (Main)"; w1; lake env lean Work/FourierU8XM/Axioms.lean 2>&1
mkdir -p /content/work/k-combine-ax
cat > /content/work/k-combine-ax/AxCheck.lean <<LEAN
import Work.GCert.Data.U8XMGlue401
import Work.GCert.Data.U8XMControls
#print axioms GXD.U8XM.exitRC
#print axioms GXD.U8XM.scalOKZ
#print axioms GXD.U8XM.validX
#print axioms GXD.U8XM.ginvX
#print axioms GXD.U8XM.ginv_cost
#print axioms GXD.U8XM.fitsS
#print axioms GXD.U8XM.price
#print axioms GXD.ginvXOfZ
#print axioms GXD.exit_of_zchecksZ
#print axioms GXD.U8XM.zcheck_control
#print axioms GXD.U8XM.classic_rejected
#print axioms GXD.U8XM.endkx_control
LEAN
echo "== axioms (glue, controls)"; w1; lake env lean /content/work/k-combine-ax/AxCheck.lean 2>&1
echo "== sorry grep"; grep -rnw sorry --include=*.lean Work/KMX Work/KLow Work/XMut Work/GCert/Data/*$N* Work/GCert/Data/Gen/$N* Work/Fourier$N Work/Fourier${N}Neg | grep -v "\`sorry\`"
echo "== forbidden constructs"; grep -rnE "^\s*axiom |native_decide|implemented_by|@\[extern|unsafe |ofReduceBool|skipKernelTC" --include=*.lean Work/KMX Work/KLow Work/XMut Work/GCert/Data/*$N* Work/GCert/Data/Gen/$N* Work/Fourier$N XExit; echo "forbidden grep done"
w0; echo "== Compare.lean positive $(date -u)"
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.Fourier$N.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N 2>&1; echo "compare rc=$?"
w0; echo "== Compare.lean NEGATIVE (must FAIL) $(date -u)"
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.Fourier${N}Neg.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N 2>&1; echo "compare-neg rc=$?"
echo "== lake build --no-build"; lake build --no-build $TG 2>&1 | tail -3
echo "CHECKS END $(date -u)"
} > $L/checks2.out 2>&1
bash $P/scripts/cmp.sh $N $T $L/cmp > $L/cmp.out 2>&1
echo "RUN2 END $(date -u)" >> $L/cmp.out
