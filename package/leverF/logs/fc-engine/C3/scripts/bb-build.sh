#!/bin/bash
# bB-data (wht45/allheads): build design B's fill-dependent modules (suffix 3945H / S5H / L3945H) in my tree on top of
# the pre-built ENG closure.  Phase A: the light kernel-data modules (core Lean, < 1.5 GB each) with LIGHT_THREADS Lean
# threads; phase B: the Mathlib modules with at most 3 heavy Lean processes on this Colab (3 minus the others, >= 1).
# Then the hypothesis-mode L3945H chain.  usage: bb-build.sh  (stage in $R/in/lean, $R/in/hyp-L3945H)
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
(cd $R/in/lean && sha256sum -c --quiet ../lean.sha256) || fail "stage hashes"
bash $S/isk-install.sh B3945H $R/in/lean || fail "install 3945H"
cd $T
L=$LOGROOT/build; mkdir -p $L
DATA="Work.XHostT.Data.S5H.Order Work.XHostT.Data.S5H.Heads $(for j in $(seq 0 129); do printf 'Work.XHostT.Data.S5H.P%d ' $j; done) $(for j in $(seq 0 26); do printf 'Work.XHostT.Data.S5H.GlobC%d ' $j; done)"
step lightA env LEAN_NUM_THREADS=${LIGHT_THREADS:-5} lake build $DATA
rc=$?; tail -2 $L/lightA.log; [ $rc = 0 ] || fail "phase A rc=$rc"
MODS="Work.XHostT.Data.S5H.Whole Work.XHostT.Data.S5H.Glob Work.XHostT.Data.S5H.DispLink Work.XHostT.Inst.S5H \
 Work.XHostT.Inst.S5HFacts Work.XHostT.Inst.S5HHB Work.Mixed.Rate3945H Work.Mixed.Rate3945HSharp Work.Mixed.R5Joint3945H \
 Work.Mixed.Bridge.Inst3945H Work.Mixed.Bridge.HeadsFirst3945H Work.Mixed.Bridge.HFTest3945HW"
wl 14000
THR=$(( MAXHEAVY - $(heavy) )); [ $THR -ge 1 ] || THR=1
echo "phase B LEAN_NUM_THREADS=$THR"
step heavyB env LEAN_NUM_THREADS=$THR lake build $MODS
rc=$?; grep -E "✖|error:" $L/heavyB.log | head -30; tail -2 $L/heavyB.log; [ $rc = 0 ] || fail "phase B rc=$rc"
# axioms of every final theorem (the modules print them; collect)
grep -hE "depends on axioms|does not depend" $L/lightA.log $L/heavyB.log | grep -E "S5H|3945H|DS5H" | sort -u > $L/axioms-3945H.txt
wc -l < $L/axioms-3945H.txt
axioms_ok $L/axioms-3945H.txt && echo "AXIOMS 3945H: standard only" || fail "axioms 3945H"
echo "BUILD 3945H OK"
