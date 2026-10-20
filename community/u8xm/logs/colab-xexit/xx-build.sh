#!/bin/bash
# c-xexit (13 Oct 2026): compile XExit.{Pack,Halves,Compat,Hub,End} inside an overlay copy of Shea's lean/complex588
# (base tree /content/work/shea, already built). Follows x-exit-lean/COLAB.sh, except: no clone/elan/cache steps
# (the base tree has them), and ONE lake invocation per module so that at most one Lean process runs at a time.
set -uo pipefail
. /content/work/env.sh
P=/content/work/c-xexit/lean/complex588
L=/content/work/logs/c-xexit
cd "$P"
waitlean() { while [ "$(pgrep -c -x lean)" -ge 2 ]; do echo "waiting: $(pgrep -c -x lean) lean processes running ($(date -u +%T))"; sleep 30; done; }
echo "== start $(date -u)"; lean --version; lake --version
echo "== pre: base closure up to date in the copy (before any change)"
lake build --no-build Work.Bridge401.End Work.Gauge.Inv Work.Bridge401.Sched; echo "PRE-NOBUILD EXIT $?"
echo "== unpack our sources"
tar xzf /content/work/xx-in/xexit-src.tar.gz && sha256sum -c XExit/SHA256SUMS; echo "SHA-OURS EXIT $?"
echo "== the five Shea files our modules import (must match the Mac's vendored copies)"
sha256sum Work/Gauge/Defs.lean Work/Gauge/Inv.lean Work/Bank401/Adapt.lean Work/Bank401/Frames.lean Work/Bridge401/Sched.lean Work/Bridge401/End.lean
grep -q 'lean_lib XExit' lakefile.lean || printf '\n-- x-exit-lean (exit-frame helpers)\nlean_lib XExit where globs := #[`XExit.+]\n' >> lakefile.lean
tail -3 lakefile.lean
echo "== post-lakefile: base closure still up to date"
lake build --no-build Work.Bridge401.End Work.Gauge.Inv Work.Bridge401.Sched; echo "POST-NOBUILD EXIT $?"
rc=0
for m in XExit.Pack XExit.Halves XExit.Compat XExit.Hub XExit.End; do
  waitlean
  echo "== build $m  $(date -u +%T)"
  LEAN_NUM_THREADS=2 /usr/bin/time -v lake build "$m" > "$L/xx-build-$m.log" 2>&1; e=$?
  cat "$L/xx-build-$m.log" | grep -vE '^\s+(Average|Socket|Signals|Swaps|File system|Page size|Voluntary|Involuntary|Minor|Major|Exit status|Percent|Average)' 
  echo "BUILD $m EXIT $e"; [ $e -eq 0 ] || rc=1
done
echo "== all five up to date"
lake build --no-build XExit.Pack XExit.Halves XExit.Compat XExit.Hub XExit.End; echo "FINAL-NOBUILD EXIT $?"
echo "== axioms"
cat > XExitAxioms.lean <<'LEAN'
import XExit.End
import XExit.Compat
import XExit.Hub
#print axioms OAI.PowerSaving.XEX.GInvX.toE'
#print axioms OAI.PowerSaving.XEX.GLabX.inv
#print axioms OAI.PowerSaving.XEX.fitsX_ofLab
#print axioms OAI.PowerSaving.XEX.HubTail.hHsplit
#print axioms OAI.PowerSaving.XEX.HubTail.hHr
#print axioms OAI.PowerSaving.XEX.wht_of_ginvx_401
#print axioms OAI.PowerSaving.XEX.engine_of_ginvx_401
#print axioms OAI.PowerSaving.XEX.engine_of_halvesX_401
#print axioms OAI.PowerSaving.XEX.wht_of_halvesX_401
#print axioms OAI.PowerSaving.GX.engine_of_ginv_401
#print axioms OAI.PowerSaving.GX.wht_of_ginv_401
#check @OAI.PowerSaving.XEX.engine_of_ginvx_401
#check @OAI.PowerSaving.XEX.wht_of_ginvx_401
LEAN
waitlean
/usr/bin/time -v lake env lean XExitAxioms.lean > "$L/xx-axioms.log" 2>&1; e=$?
cat "$L/xx-axioms.log"; echo "AXIOMS EXIT $e"
echo "== negative control: the same sorryAx test must fire on a module that uses sorry"
printf 'import XExit.Pack\ntheorem xxNegCtrl : (1:Nat) = 2 := sorry\n#print axioms xxNegCtrl\n' > XExitNeg.lean
waitlean
lake env lean XExitNeg.lean > "$L/xx-neg.log" 2>&1; echo "NEG EXIT $?" >> "$L/xx-neg.log"; cat "$L/xx-neg.log"
echo "== sorry/escape grep in our sources"
grep -nwE 'sorry|admit|axiom|native_decide|implemented_by|extern|unsafe|ofReduceBool|skipKernelTC|decide_by_eval|debug.skipKernelTC' XExit/*.lean | grep -vE '^\S+:[0-9]+:\s*(--|/-|#print axioms)' || echo "no hits"
grep -n 'sorry' XExit/*.lean || echo "no 'sorry' substring at all"
echo "sorryAx in axioms log: $(grep -c sorryAx "$L/xx-axioms.log")"
echo "'declaration uses sorry' in build logs: $(cat $L/xx-build-*.log | grep -c 'declaration uses .sorry.')"
echo "== end $(date -u)  OVERALL rc=$rc"
