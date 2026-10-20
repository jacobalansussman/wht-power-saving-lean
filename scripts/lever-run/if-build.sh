#!/bin/bash
# i-final (wht45/lever): build the whole L3945 closure in my tree /content/work/i-final/T before the fixed final run.
# Phase A: the light kernel-data modules (core Lean, < 1 GB each) with 8 Lean threads; phase B: everything else with
# at most 3 heavy Lean processes on this Colab (LEAN_NUM_THREADS = 3 minus the others).
set -uo pipefail
export ISK_R=/content/work/i-final ISK_LOG=/content/work/logs/i-final ISK_SCRIPTS=/content/work/i-final/scripts
. $ISK_SCRIPTS/isk-lib.sh
G=$R/in/final-L3945-run
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "stage hashes"
(cd $R/in && sha256sum -c --quiet ENG.sha256.chk) || fail "ENG hashes"
bash $S/isk-tree.sh || fail tree
bash $S/isk-install.sh L3945-extra $R/in/engine || fail "install extra"
bash $S/isk-install.sh L3945 $G/stage || fail install
cd $T
L=$LOGROOT/build-a; mkdir -p $L
step lightA env LEAN_NUM_THREADS=8 lake build Work.XHostT.Data.S5.P0 Work.XHostT.Data.S5.P1 Work.XHostT.Data.S5.P2 Work.XHostT.Data.S5.P3 Work.XHostT.Data.S5.P4 Work.XHostT.Data.S5.P5 Work.XHostT.Data.S5.P6 Work.XHostT.Data.S5.P7 Work.XHostT.Data.S5.P8 Work.XHostT.Data.S5.P9 Work.XHostT.Data.S5.P10 Work.XHostT.Data.S5.P11 Work.XHostT.Data.S5.P12 Work.XHostT.Data.S5.P13 Work.XHostT.Data.S5.P14 Work.XHostT.Data.S5.P15 Work.XHostT.Data.S5.P16 Work.XHostT.Data.S5.P17 Work.XHostT.Data.S5.P18 Work.XHostT.Data.S5.P19 Work.XHostT.Data.S5.P20 Work.XHostT.Data.S5.P21 Work.XHostT.Data.S5.P22 Work.XHostT.Data.S5.P23 Work.XHostT.Data.S5.P24 Work.XHostT.Data.S5.P25 Work.XHostT.Data.S5.P26 Work.XHostT.Data.S5.P27 Work.XHostT.Data.S5.P28 Work.XHostT.Data.S5.P29 Work.XHostT.Data.S5.P30 Work.XHostT.Data.S5.P31 Work.XHostT.Data.S5.P32 Work.XHostT.Data.S5.P33 Work.XHostT.Data.S5.P34 Work.XHostT.Data.S5.P35 Work.XHostT.Data.S5.P36 Work.XHostT.Data.S5.P37 Work.XHostT.Data.S5.P38 Work.XHostT.Data.S5.P39 Work.XHostT.Data.S5.P40 Work.XHostT.Data.S5.P41 Work.XHostT.Data.S5.P42 Work.XHostT.Data.S5.P43 Work.XHostT.Data.S5.P44 Work.XHostT.Data.S5.P45 Work.XHostT.Data.S5.P46 Work.XHostT.Data.S5.P47 Work.XHostT.Data.S5.P48 Work.XHostT.Data.S5.P49 Work.XHostT.Data.S5.P50 Work.XHostT.Data.S5.P51 Work.XHostT.Data.S5.P52 Work.XHostT.Data.S5.P53 Work.XHostT.Data.S5.P54 Work.XHostT.Data.S5.P55 Work.XHostT.Data.S5.P56 Work.XHostT.Data.S5.P57 Work.XHostT.Data.S5.P58 Work.XHostT.Data.S5.P59 Work.XHostT.Data.S5.P60 Work.XHostT.Data.S5.P61 Work.XHostT.Data.S5.P62 Work.XHostT.Data.S5.P63 Work.XHostT.Data.S5.P64 Work.XHostT.Data.S5.P65 Work.XHostT.Data.S5.P66 Work.XHostT.Data.S5.P67 Work.XHostT.Data.S5.P68 Work.XHostT.Data.S5.P69 Work.XHostT.Data.S5.P70 Work.XHostT.Data.S5.P71 Work.XHostT.Data.S5.P72 Work.XHostT.Data.S5.P73 Work.XHostT.Data.S5.P74 Work.XHostT.Data.S5.P75 Work.XHostT.Data.S5.P76 Work.XHostT.Data.S5.P77 Work.XHostT.Data.S5.P78 Work.XHostT.Data.S5.P79 Work.XHostT.Data.S5.P80 Work.XHostT.Data.S5.P81 Work.XHostT.Data.S5.P82 Work.XHostT.Data.S5.P83 Work.XHostT.Data.S5.P84 Work.XHostT.Data.S5.P85 Work.XHostT.Data.S5.P86 Work.XHostT.Data.S5.P87 Work.XHostT.Data.S5.P88 Work.XHostT.Data.S5.P89 Work.XHostT.Data.S5.P90 Work.XHostT.Data.S5.P91 Work.XHostT.Data.S5.P92 Work.XHostT.Data.S5.P93 Work.XHostT.Data.S5.P94 Work.XHostT.Data.S5.P95 Work.XHostT.Data.S5.P96 Work.XHostT.Data.S5.P97 Work.XHostT.Data.S5.P98 Work.XHostT.Data.S5.P99 Work.XHostT.Data.S5.P100 Work.XHostT.Data.S5.P101 Work.XHostT.Data.S5.P102 Work.XHostT.Data.S5.P103 Work.XHostT.Data.S5.P104 Work.XHostT.Data.S5.P105 Work.XHostT.Data.S5.P106 Work.XHostT.Data.S5.P107 Work.XHostT.Data.S5.P108 Work.XHostT.Data.S5.P109 Work.XHostT.Data.S5.P110 Work.XHostT.Data.S5.P111 Work.XHostT.Data.S5.P112 Work.XHostT.Data.S5.P113 Work.XHostT.Data.S5.P114 Work.XHostT.Data.S5.P115 Work.XHostT.Data.S5.P116 Work.XHostT.Data.S5.P117 Work.XHostT.Data.S5.P118 Work.XHostT.Data.S5.P119 Work.XHostT.Data.S5.P120 Work.XHostT.Data.S5.P121 Work.XHostT.Data.S5.P122 Work.XHostT.Data.S5.P123 Work.XHostT.Data.S5.P124 Work.XHostT.Data.S5.P125 Work.XHostT.Data.S5.P126 Work.XHostT.Data.S5.P127 Work.XHostT.Data.S5.P128 Work.XHostT.Data.S5.P129 Work.XHostT.Data.S5.GlobC0 Work.XHostT.Data.S5.GlobC1 Work.XHostT.Data.S5.GlobC2 Work.XHostT.Data.S5.GlobC3 Work.XHostT.Data.S5.GlobC4 Work.XHostT.Data.S5.GlobC5 Work.XHostT.Data.S5.GlobC6 Work.XHostT.Data.S5.GlobC7 Work.XHostT.Data.S5.GlobC8 Work.XHostT.Data.S5.GlobC9 Work.XHostT.Data.S5.GlobC10 Work.XHostT.Data.S5.GlobC11 Work.XHostT.Data.S5.GlobC12 Work.XHostT.Data.S5.GlobC13 Work.XHostT.Data.S5.GlobC14 Work.XHostT.Data.S5.GlobC15 Work.XHostT.Data.S5.GlobC16 Work.XHostT.Data.S5.GlobC17 Work.XHostT.Data.S5.GlobC18 Work.XHostT.Data.S5.GlobC19 Work.XHostT.Data.S5.GlobC20 Work.XHostT.Data.S5.GlobC21 Work.XHostT.Data.S5.GlobC22 Work.XHostT.Data.S5.GlobC23 Work.XHostT.Data.S5.GlobC24 Work.XHostT.Data.S5.GlobC25 Work.XHostT.Data.S5.GlobC26 Work.XHostT.Data.S5.Order Work.XHostT.Data.S5.Heads
tail -3 $L/lightA.log
L=$LOGROOT/build-b; mkdir -p $L
wl 20000
THR=$(( MAXHEAVY - $(heavy) )); [ $THR -ge 1 ] || THR=1
echo "phase B LEAN_NUM_THREADS=$THR"
step engine env LEAN_NUM_THREADS=$THR lake build Work.L3945.Engine $(cat $G/meta/targets.txt)
rc=$?
grep -E "✖|error" $L/engine.log | head -40
tail -3 $L/engine.log
echo "PHASE B rc=$rc"
