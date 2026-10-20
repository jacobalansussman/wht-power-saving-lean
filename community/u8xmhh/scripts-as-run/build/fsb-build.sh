#!/bin/bash
# fs-build (wht43 bound-check, 13 Oct 2026, Colab #2): build and check GXD.U8XMHH in MY tree (fsb-setup.sh made it).
# Lake route as root (my copy's .lake is root-owned after setup; traces written so the comparator has nothing to
# rebuild), ONE Lake job at a time (LEAN_NUM_THREADS=1), each Lean step only when fewer than 2 Lean processes over 2 GB
# run on the machine and enough memory is free.  Then: axioms, sorry/forbidden grep, challenge derivation, hostSchedK
# negative, Compare.lean positive and negative.  The official comparator is fsb-cmp.sh (after the chown to 65534).
# usage: fsb-build.sh      logs /content/work/logs/fs-build/
set -uo pipefail
. /content/work/env.sh
F=/content/work/fs-build; T=$F/B/lean/complex588; L=/content/work/logs/fs-build; mkdir -p $L $F/ax
N=U8XMHH; POS=16197705; NEG=16197706
REF=/content/work/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean
[ -e $F/SETUP-OK ] || { echo "no SETUP-OK"; exit 1; }
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
avail() { free -m | awk '/Mem:/ {print $7}'; }
wl() { local need=${1:-12000} w=0; while [ "$(heavy)" -ge 2 ] || [ "$(avail)" -lt $need ]; do sleep 10; w=$((w+10)); done; echo "   (slot: waited ${w}s; heavy lean now $(heavy), avail $(avail) MB)"; }
step() { local lab=$1; shift; local s=$(date +%s); echo "=== $lab START $(date -u +%T)"
  /usr/bin/time -v -o $L/$lab.time "$@" > $L/$lab.log 2>&1; local rc=$?
  echo "=== $lab END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$lab.time) $(date -u +%T)"; return $rc; }
lk() { local lab=$1; shift; wl 12000; step $lab env LEAN_NUM_THREADS=1 lake build "$@"; local rc=$?
  grep -E "✖|error:|sorry|Build completed" $L/$lab.log | head -20
  grep -E "Built .* \([0-9.]+(s|m)\)" $L/$lab.log | sed -E 's/.*Built ([^ ]+) \(([0-9.]+)(s|m)\).*/   \2\3 \1/' | sort -g -r | head -8
  [ $rc -eq 0 ] || { echo "BUILD FAILED at $lab"; exit 14; }; }
cd $T
echo "FSB-BUILD START $(date -u) exponent 1 - $POS/10^10 (kernel 1 - $NEG/10^10; negative challenge 1 - $NEG/10^10)"
echo "owner of .lake: $(stat -c %U .lake) packages: $(stat -c %U .lake/packages/mathlib)"
[ "$(stat -c %U .lake/packages/mathlib)" = root ] || { echo "packages not root-owned: refusing root lake"; exit 2; }
## 1. staged Lake builds (fail fast: rate fact, host table check, network, glue, Fourier chain)
lk b1-rate Work.XHost.Data.RateU8XMHHx
lk b2-hostk Work.XHost.HSchedK Work.XHost.HPriceK
lk b3-hosttab Work.XHost.Data.HU8XMSched Work.XHost.Data.InstU8XMHH
lk b4-network Work.XHost.HBank Work.XHost.HEnd
lk b5-glue Work.XHost.U8XMHH
TG="Work.FourierU8XMHH.Main Work.FourierU8XMHH.UniformFourierChallenge Work.FourierU8XMHHNeg.UniformFourierChallenge Work.FourierU8XMHH.Axioms"
lk b6-fourier $TG
wl 8000; step nobuild lake build --no-build $TG Work.FourierU8XM.Main Work.XHost.U8XMHH; tail -1 $L/nobuild.log
echo "'declaration uses sorry' in the builds: $(cat $L/b?-*.log | grep -c "declaration uses 'sorry'") (the two in the challenge files are the challenge statements)"
cat $L/b?-*.log | grep "declaration uses 'sorry'"
## 2. axioms
wl 12000; step axioms lake env lean Work/Fourier$N/Axioms.lean; cat $L/axioms.log
cat > $F/ax/AxCheckU8XMHH.lean <<'LEAN'
import Work.FourierU8XMHH.Seam
#print axioms GXD.U8XMHH.hs_ok
#print axioms GXD.U8XMHH.hfits
#print axioms GXD.U8XMHH.hwidth
#print axioms GXD.U8XMHH.hprice
#print axioms OAI.PowerSaving.BridgeRate.U8XMHH.unit_cost_x
#print axioms OAI.PowerSaving.BlockAccounting.fold_U8XMHHx
#print axioms OAI.PowerSaving.XHost.HostSched.ofTab_fitsR
#print axioms OAI.PowerSaving.XHost.engine_of_ginv_host
#print axioms OAI.PowerSaving.XHost.wht_of_ginv_host
#print axioms OAI.PowerSaving.SeamU8XMHH.engine_U8XMHHx
#print axioms OAI.PowerSaving.RAM.hillsU8XMHH_program
#print axioms GXD.U8XM.ginv_cost
#print axioms GXD.U8XM.hent
#print axioms GXD.U8XM.exitRC
#print axioms GXD.U8XM.scalOKZ
#print axioms GXD.U8XM.validX
#print axioms GXD.ginvXOfZ
#print axioms GXD.exit_of_zchecksZ
#check @OAI.PowerSaving.XHost.engine_of_ginv_host
#check @OAI.PowerSaving.SeamU8XMHH.engine_U8XMHHx
LEAN
wl 12000; step axioms-glue lake env lean $F/ax/AxCheckU8XMHH.lean; cat $L/axioms-glue.log
## 3. hostSchedK negative control (fs-price's one-off copy of the table: one remainder off by one -> false)
cp $F/in/HsU8XMHHNeg.lean $F/ax/HsU8XMHHNeg.lean
wl 8000; step hsneg lake env lean $F/ax/HsU8XMHHNeg.lean; echo "hostSchedK negative rc=$? (0 = the kernel proved '= false')"; tail -3 $L/hsneg.log
## 4. source checks: every file I added (XHost network, U8XMHH data/glue, Fourier chain) and the U8XM overlay
(cd $T && { find Work/XHost Work/FourierU8XMHH Work/FourierU8XMHHNeg -name '*.lean'; sed -n 's#^[0-9a-f]*  overlay/##p' /content/work/kb/pkg/B/k-combine-u8xm/SHA256SUMS | grep '\.lean$'; } | sort -u) > $L/checked-lean-files.txt
echo "files checked: $(wc -l < $L/checked-lean-files.txt)"
(cd $T && xargs grep -nw sorry < $L/checked-lean-files.txt | grep -v '`sorry`') > $L/sorry-grep.txt; echo "sorry lines:"; cat $L/sorry-grep.txt
(cd $T && xargs grep -nwE 'admit' < $L/checked-lean-files.txt) > $L/admit-grep.txt; echo "admit lines: $(wc -l < $L/admit-grep.txt)"
(cd $T && xargs grep -nE '^\s*(private |protected |noncomputable )?axiom |native_decide|implemented_by|@\[extern|unsafe |ofReduceBool|ofReduceNat|skipKernelTC|debug\.skipKernelTC|Lean\.Elab|elab_rules|macro_rules' < $L/checked-lean-files.txt) > $L/forbidden-grep.txt; echo "forbidden-construct lines: $(wc -l < $L/forbidden-grep.txt)"; head -20 $L/forbidden-grep.txt
## 5. challenge statements re-derived from OpenAI's original
echo "OpenAI ref: $(sha256sum $REF)"
for spec in "$POS Fourier$N" "$NEG Fourier${N}Neg"; do set -- $spec
  bash $F/scripts/derive_challenges.sh $REF $1 $N > $L/derived-$2.lean
  cmp -s $L/derived-$2.lean $T/Work/$2/UniformFourierChallenge.lean && echo "challenge $2 (1 - $1/10^10): byte-identical to my derivation" || echo "CHALLENGE $2 DIFFERS FROM MY DERIVATION"; done
echo "Compare.lean: $(sha256sum tools/Compare.lean) (ours: $(sha256sum /content/work/ours/Compare.lean | cut -c1-64))"
## 6. Compare.lean positive, then negative (must FAIL)
wl 20000; step compare-pos lake env lean --run tools/Compare.lean Work.Fourier$N.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N; echo "compare-pos rc=$?"; grep -E "RESULT" $L/compare-pos.log
wl 20000; step compare-neg lake env lean --run tools/Compare.lean Work.Fourier${N}Neg.UniformFourierChallenge Work.Fourier$N.Main OAI.PowerSaving.transform_main$N OAI.PowerSaving.convolution_main$N; echo "compare-neg rc=$? (must be nonzero)"; grep -E "RESULT|differs" $L/compare-neg.log | head -5
touch $F/BUILD-CHECKS-DONE
echo "FSB-BUILD END $(date -u)"
