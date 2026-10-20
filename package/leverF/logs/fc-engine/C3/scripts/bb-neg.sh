#!/bin/bash
# bB-data: NEGATIVE controls of design B's data in Lean.  Each NegSome / Neg* file must FAIL (the kernel refutes the
# claim), each NegNone file must compile (the checker itself returns none on the perturbed table).
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
bash $S/isk-install.sh NEG $R/in/neg || fail "install neg"
L=$LOGROOT/neg; mkdir -p $L; cd $T
step negcalls env LEAN_NUM_THREADS=2 lake build Work.XHostT.Data.S5Hn.Calls Work.XHostT.Data.S5Hw.Calls || fail "neg Calls"
ok=1
for f in Work/XHostT/Data/S5Hn/NegSome*.lean Work/XHostT/Data/S5Hw/NegSome*.lean Work/XHostT/Data/S5Hneg/NegHeads.lean \
         Work/XHostT/Data/S5Hneg/NegHeadsOld.lean Work/XHostT/Data/S5Hneg/NegRes.lean Work/XHostT/Data/S5Hneg/NegGlob.lean; do
  b=$(echo $f | tr / _); step neg-$b lake env lean $f; rc=$?
  if [ $rc != 0 ] && grep -qE "error" $L/neg-$b.log; then echo "NEG-FAILS-AS-REQUIRED $f: $(grep -m1 -oE 'decide failed.{0,80}|proposition.{0,60}false|error.{0,100}' $L/neg-$b.log)"; else echo "NEGATIVE DID NOT FAIL: $f"; ok=0; fi
done
for f in Work/XHostT/Data/S5Hn/NegNone*.lean Work/XHostT/Data/S5Hw/NegNone*.lean; do
  b=$(echo $f | tr / _); step neg-$b lake env lean $f; rc=$?
  if [ $rc = 0 ] && ! grep -q "error" $L/neg-$b.log; then echo "NONE-COMPILES $f"; else echo "NegNone did not compile: $f"; ok=0; fi
done
[ $ok = 1 ] && echo "NEG ALL AS REQUIRED" || fail "negatives"
