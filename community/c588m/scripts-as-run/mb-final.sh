#!/bin/bash
# m-build (13 Oct 2026): the checks after the 70 C588M modules are built (replaces the tail of build-c588m-mb.sh /
# mb-build.sh, which released the Lean slot between steps and starved on a crowded machine). Waits ONCE for a free
# Lean slot (3 consecutive samples of `pgrep -c lean` < 2), then runs its heavy steps back to back, one at a time:
#   axiom report; Compare.lean (ours, byte-identical to the tree's tools/Compare.lean) positive and negative;
#   the official comparator (c-base's cb-cmp.sh procedure in jail2.sh) positive and negative.
. /content/work/env.sh
T=/content/work/mb/shea/lean/complex588; L=/content/work/logs/m-build; cd $T || exit 90
export LEAN_NUM_THREADS=2
step() { # step <name> <command...>
  local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"
  return $rc
}
waitslot() { local ok=0; while [ $ok -lt 3 ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 10; else ok=$((ok+1)); sleep 3; fi; done; }
T0=$(date +%s)
echo "FINAL START $(date -u)"
/content/work/cb-sampler.sh $L/mem-sampler-final.log & SPID=$!
trap 'kill $SPID 2>/dev/null' EXIT
echo "Compare.lean: $(sha256sum /content/work/ours/Compare.lean tools/Compare.lean | tr '\n' ' ')"
if [ "${HANDOFF:-0}" = 1 ]; then
  # hand-over: mb-build.sh (tmux mbbuild) holds our Lean slot and is running axioms/compare; take the slot over
  # from it when its checks are done (no release in between), then stop it before its redundant steps.
  until grep -qE "ALL CHECKS PASSED|SOME CHECK FAILED|BUILD END" $L/mb-build.out; do sleep 3; done
  tmux -L mbbuild kill-server 2>/dev/null; sleep 2
  echo "HANDOFF from mb-build.sh $(date -u); its result: $(grep -E 'ALL CHECKS PASSED|SOME CHECK FAILED' $L/mb-build.out)"
  grep -E "^(PASS|FAIL) " $L/mb-build.out
else
waitslot; echo "SLOT ACQUIRED $(date -u) after $(( $(date +%s)-T0 ))s; lean processes now: $(pgrep -c lean)"
fi
if [ "${FROM:-all}" != cmp ]; then
step axioms lake env lean Work/FourierC588M/Axioms.lean
step compare-fourier lake env lean --run /content/work/ours/Compare.lean Work.FourierC588M.UniformFourierChallenge \
  Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
step compare-negative lake env lean --run /content/work/ours/Compare.lean Work.FourierC588MNeg.UniformFourierChallenge \
  Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
NEG=$?
echo "(negative control must be nonzero: rc=$NEG)"
python3 -I - "$L" "$NEG" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1]); neg = int(sys.argv[2])
ok = []
ok.append(('positive comparison RESULT: PASS', 'RESULT: PASS' in (p/'compare-fourier.log').read_text()))
bad = (p/'compare-negative.log').read_text()
ok.append(('negative comparison fails (rc != 0)', neg != 0))
ok.append(('negative: definition body differs, decimalExponentC588M', 'definition body differs' in bad and 'decimalExponentC588M' in bad))
ax = (p/'axioms.log').read_text()
ok.append(('no sorryAx', 'sorryAx' not in ax))
ok.append(('three theorems on [propext, Classical.choice, Quot.sound]', ax.count('[propext, Classical.choice, Quot.sound]') == 3))
for name, v in ok: print(('PASS ' if v else 'FAIL ') + name)
print('ALL CHECKS PASSED' if all(v for _, v in ok) else 'SOME CHECK FAILED')
PY
fi
chown -R 65534:65534 $T/.lake
bash /content/work/mb/mb-cmp.sh m-build-A1 $T comparator/fourier-c588m.json > $L/cmp-A1.out 2>&1; echo "COMPARATOR A1 positive rc=$? $(date -u)"
grep -E "accepts the solution|Your solution is okay|uncaught exception|=== .* END" $L/cmp-A1.out
bash /content/work/mb/mb-cmp.sh m-build-A2neg $T comparator/fourier-c588m-neg.json > $L/cmp-A2neg.out 2>&1; echo "COMPARATOR A2 negative rc=$? (must be nonzero) $(date -u)"
grep -E "accepts the solution|Your solution is okay|uncaught exception|=== .* END" $L/cmp-A2neg.out
step nobuild lake build --no-build Work.FourierC588M.Main Work.FourierC588M.UniformFourierChallenge \
  Work.FourierC588MNeg.UniformFourierChallenge Work.GCert.Data.C588MControls Work.FourierC588M.Axioms
PKG=/content/work/m-gen/c588m
{
  echo "# grep of the overlay files in the tree for sorry/admit/axiom/native_decide/implemented_by/extern/unsafe/opaque (word match, comments included)"
  (cd $PKG/overlay && find . -name '*.lean' | sed 's|^\./||' | sort) | while read -r f; do
    grep -nwE 'sorry|admit|axiom|native_decide|implemented_by|extern|unsafe|opaque|sorryAx' "$T/$f" | sed "s|^|$f:|"
  done
  echo "# declarations starting with 'axiom' anywhere under OAI Work WHTCheck:"
  grep -rnE '^\s*(private |protected )?axiom ' OAI Work WHTCheck | head -20
  echo "# tracked files of the tree changed vs git HEAD:"
  git status --short --untracked-files=no
  echo "# untracked files: $(git status --short --untracked-files=all | grep -c '^??')"
} > $L/source-scan.txt 2>&1
echo "FINAL END $(date -u) total_wall=$(( $(date +%s)-T0 ))s"
awk '{ if ($5+0>m) m=$5+0; if ($4+0>s) s=$4+0 } END { print "sampler: max single process RSS " m " MB, max summed lean/lake RSS " s " MB (all agents)" }' $L/mem-sampler-final.log
