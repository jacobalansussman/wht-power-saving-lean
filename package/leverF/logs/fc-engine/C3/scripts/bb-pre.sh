#!/bin/bash
# bB-data (wht45/allheads): my tree on this Colab = copy of ours/T + i-final's ENG overlay (fill-independent engine,
# network, bridge modules); then build the fill-independent modules my 3945H files import.  Lean threads = 3 minus
# the other heavy Lean processes (at least 1); waits for 14 GB free.
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
(cd $R/in && sha256sum -c --quiet ENG.sha256.chk) || fail "ENG hashes"
bash $S/isk-tree.sh || fail tree
bash $S/isk-install.sh ENG $R/in/engine || fail "install ENG"
cd $T
L=$LOGROOT/pre; mkdir -p $L
wl 14000
THR=$(( MAXHEAVY - $(heavy) )); [ $THR -ge 1 ] || THR=1
echo "pre LEAN_NUM_THREADS=$THR"
step pre env LEAN_NUM_THREADS=$THR lake build Work.XHostT.HFW2979.GeomInstW Work.Mixed.Bridge.HFTest Work.Mixed.Bridge.Inst \
  Work.Mixed.Bridge.HeadsFirst Work.XHostT.TTab Work.XHostT.TEnd Work.XHostT.KHead Work.XHostT.KHeadC Work.XHostT.KHist \
  Work.XHostT.KHistC Work.XHostT.KWhole Work.XHostT.KGlue Work.XHostT.Data.D3945.Disp Work.Mixed.TwoTarget \
  Work.BlockAccounting.Sharp Work.BlockAccounting.Moment4 Work.Mixed.SeamConst Work.XHostT.HF.GeomInst Work.Mixed.HeadM1
rc=$?
grep -E "✖|error" $L/pre.log | head -20; tail -3 $L/pre.log
echo "PRE rc=$rc"
