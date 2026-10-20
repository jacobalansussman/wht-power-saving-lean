#!/bin/bash
# v-fs-rebuild (skeptic, Colab #1): build the 32 new modules (12 XHost network + 20 U8XMHH) in MY copy of k-combine's
# original U8XM tree, ONE module per Lake call in topological order (so at most one new Lean process at a time),
# LEAN_NUM_THREADS=1; then axioms, usage closure, Compare.lean positive and negative.
set -uo pipefail
. /content/work/env.sh
D=/content/work/v-fs-rebuild; T=$D/T; L=/content/work/logs/v-fs-rebuild
[ -e $D/SETUP-OK ] || { echo "no SETUP-OK"; exit 1; }
cd $T
[ "$(stat -c %U .lake/packages/mathlib)" = root ] || { echo "packages not root-owned: refusing root lake"; exit 2; }
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
wl() { local w=0; while [ "$(heavy)" -ge 2 ] || [ "$(free -m | awk '/Mem:/ {print $7}')" -lt 14000 ]; do sleep 10; w=$((w+10)); done; [ $w -gt 0 ] && echo "   (waited ${w}s for a slot)"; }
step() { local lab=$1; shift; local s=$(date +%s)
  /usr/bin/time -v -o $L/$lab.time "$@" > $L/$lab.log 2>&1; local rc=$?
  echo "=== $lab rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$lab.time) $(date -u +%T)"; return $rc; }
echo "VFR-BUILD START $(date -u)"
echo "olean before: XHost $(find .lake/build/lib/lean/Work/XHost -name '*.olean' 2>/dev/null | wc -l), FourierU8XMHH $(find .lake/build/lib/lean/Work/FourierU8XMHH -name '*.olean' 2>/dev/null | wc -l)"
MODS="Work.XHost.HostK Work.XHost.HRoles Work.XHost.HDef Work.XHost.HStage Work.XHost.HCert Work.XHost.Frames Work.XHost.HSched Work.XHost.HGeom Work.XHost.HEnd Work.XHost.HBank Work.XHost.HSchedK Work.XHost.HPriceK Work.XHost.Data.HU8XMSched Work.XHost.Data.RateU8XMHHx Work.XHost.Data.InstU8XMHH Work.XHost.U8XMHH Work.FourierU8XMHH.Seam Work.FourierU8XMHH.SectorAlgorithm Work.FourierU8XMHH.SynchronizedAlgorithm Work.FourierU8XMHH.WorkingTransform Work.FourierU8XMHH.WorkingCompiler Work.FourierU8XMHH.WorkingPreparation Work.FourierU8XMHH.ArbitraryLength Work.FourierU8XMHH.TransformProgram Work.FourierU8XMHH.ConvolutionProgram Work.FourierU8XMHH.UniformBounds Work.FourierU8XMHH.Asymptotics Work.FourierU8XMHH.Goal Work.FourierU8XMHH.Main Work.FourierU8XMHH.Axioms Work.FourierU8XMHH.UniformFourierChallenge Work.FourierU8XMHHNeg.UniformFourierChallenge"
i=0
for M in $MODS; do i=$((i+1)); wl
  step "m$(printf %02d $i)-$M" env LEAN_NUM_THREADS=1 lake build $M || { echo "BUILD FAILED at $M"; tail -30 $L/m$(printf %02d $i)-$M.log; exit 14; }
  grep -E "^(✔|⚠|✖)|Built|Replayed" $L/m$(printf %02d $i)-$M.log | grep -v "Replayed" | cut -c1-160 | tail -4
done
step nobuild lake build --no-build $MODS Work.FourierU8XM.Main; echo "no-build: $(tail -1 $L/nobuild.log)"
echo "'declaration uses' sorry warnings in my builds:"; cat $L/m*.log | grep -n "declaration uses 'sorry'" | cut -c1-200
cat $L/m*.log | grep -E "^✖|error:" | head
echo "olean after: XHost $(find .lake/build/lib/lean/Work/XHost -name '*.olean' | wc -l), FourierU8XMHH $(find .lake/build/lib/lean/Work/FourierU8XMHH -name '*.olean' | wc -l)"
## axioms (Axioms.lean prints during build; also directly)
cat $L/m30-Work.FourierU8XMHH.Axioms.log | grep -A3 "depends on" | head -20
## Compare.lean (our stand-in), positive then negative
echo "Compare.lean $(sha256sum tools/Compare.lean)"
N=U8XMHH
wl; step compare-pos lake env lean --run tools/Compare.lean Work.Fourier$N.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N; echo "compare-pos rc=$?"; grep -E "RESULT|FAIL|OK:|axioms:" $L/compare-pos.log | cut -c1-250
wl; step compare-neg lake env lean --run tools/Compare.lean Work.Fourier${N}Neg.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N; echo "compare-neg rc=$? (must be nonzero)"; grep -E "RESULT|FAIL" $L/compare-neg.log | head -6 | cut -c1-250
## usage closure
wl; step uses lake env lean --run $D/scripts/Uses.lean Work.FourierU8XMHH.Main OAI.PowerSaving.transform_mainU8XMHH OAI.PowerSaving.convolution_mainU8XMHH; echo "uses rc=$?"; cat $L/uses.log | cut -c1-400
echo "VFR-BUILD END $(date -u)"
