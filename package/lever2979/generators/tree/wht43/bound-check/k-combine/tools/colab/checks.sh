#!/bin/bash
# k-combine (k-lower's checks.sh for the stacked route): axioms of the three theorems and of the glue facts, sorry
# grep, the two challenges, Compare.lean positive and one-unit-stronger negative, the Lean negative controls.
# usage: checks.sh <NAME> <tree> <log dir>
set -uo pipefail
. /content/work/env.sh
N=$1; T=$2; L=$3
mkdir -p $L
cd $T
o() { echo .lake/build/lib/lean/$(echo $1 | tr . /); }
cc() { while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done; f=$(echo $1 | tr . /).lean; mkdir -p $(dirname $(o $1));
       LEAN_NUM_THREADS=2 lake env lean -DautoImplicit=false -R . -o $(o $1).olean -i $(o $1).ilean $f; }
w1() { while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 2 ]; do sleep 5; done; }
echo "== axioms (Main)"; w1; LEAN_NUM_THREADS=2 lake env lean -DautoImplicit=false Work/Fourier$N/Axioms.lean 2>&1 | tee $L/axioms.log
mkdir -p /tmp/k-combine-ax
cat > /tmp/k-combine-ax/AxCheck$N.lean <<LEAN
import Work.GCert.Data.${N}Glue401
#print axioms GXD.$N.exitRC
#print axioms GXD.$N.scalOKZ
#print axioms GXD.$N.validX
#print axioms GXD.$N.ginvX
#print axioms GXD.$N.ginv_cost
#print axioms GXD.$N.fitsS
#print axioms GXD.$N.price
#print axioms GXD.ginvXOfZ
#print axioms GXD.exit_of_zchecksZ
LEAN
echo "== axioms (glue)"; w1; LEAN_NUM_THREADS=2 lake env lean -DautoImplicit=false /tmp/k-combine-ax/AxCheck$N.lean 2>&1 | tee $L/axioms-glue.log
echo "== sorry grep"; grep -rnw sorry --include=*.lean Work/KMX Work/KLow Work/XMut Work/GCert/Data/*$N* Work/GCert/Data/Gen/$N* Work/Fourier$N Work/Fourier${N}Neg | grep -v '`sorry`' | tee $L/sorry-grep.log
echo "== forbidden constructs"; grep -rnE "^\s*axiom |native_decide|implemented_by|@\[extern|unsafe |ofReduceBool|skipKernelTC" --include=*.lean Work/KMX Work/KLow Work/XMut Work/GCert/Data/*$N* Work/GCert/Data/Gen/$N* Work/Fourier$N | tee $L/forbidden-grep.log; echo "forbidden lines: $(wc -l < $L/forbidden-grep.log)"
echo "== challenges"; cc Work.Fourier$N.UniformFourierChallenge > $L/challenge.log 2>&1; echo "challenge rc=$?"
cc Work.Fourier${N}Neg.UniformFourierChallenge > $L/challenge-neg.log 2>&1; echo "challenge-neg rc=$?"
if [ -f Work/GCert/Data/${N}Controls.lean ]; then cc Work.GCert.Data.${N}Controls > $L/controls.log 2>&1; echo "controls rc=$?"; grep "depends on axioms" $L/controls.log; fi
while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 1 ]; do sleep 5; done
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.Fourier$N.UniformFourierChallenge Work.Fourier$N.Main \
  OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N > $L/compare.log 2>&1; echo "compare rc=$?"; grep RESULT $L/compare.log
while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 1 ]; do sleep 5; done
/usr/bin/time -v lake env lean --run tools/Compare.lean Work.Fourier${N}Neg.UniformFourierChallenge Work.Fourier$N.Main \
  OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N > $L/compare-neg.log 2>&1; echo "compare-neg rc=$? (must FAIL)"; grep RESULT $L/compare-neg.log
echo "CHECKS END $(date -u)"
