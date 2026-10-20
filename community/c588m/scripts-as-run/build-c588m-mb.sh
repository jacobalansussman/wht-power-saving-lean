#!/bin/bash
# build-c588m-mb.sh [tree] [label]  (m-build's copy of m-gen's build-c588m.sh: fixed PKG path, robust Lean-slot wait)
# build-c588m.sh [tree] [label]   (run inside a fresh tmux server: tmux -L m-c588m new -d -s b '...')
# m-gen (13 Oct 2026).  Builds the 78 new modules of the copy-free #588 certificate C588M ON TOP of Ryan Shea's built
# lean/complex588 tree (any variant applied), ONE module per `lake build` in the topological order of MODULES.txt (so at most one of our
# Lean processes runs at a time), waiting while 2 or more Lean processes run on the machine; then the axiom report
# and the two statement comparisons (Compare.lean): positive at 13552606/10^10, negative one unit stronger.
# Every step wrapped in /usr/bin/time -v (as c-base's cb-build.sh); logs in /content/work/logs/<label>/.
. /content/work/env.sh
P=${1:-/content/work/shea/lean/complex588}; LBL=${2:-m-c588m}
PKG=${PKG:-/content/work/m-gen/c588m}   # m-build copy: package path fixed (this copy lives in /content/work/mb)
L=/content/work/logs/$LBL; mkdir -p $L
cd $P || exit 90
export LEAN_NUM_THREADS=${LEAN_NUM_THREADS:-2}
step() { # step <name> <command...>
  local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"
  return $rc
}
waitlean() { local ok=0; while [ $ok -lt 3 ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 15; else ok=$((ok+1)); sleep 4; fi; done; }   # m-build: 3 consecutive samples below 2 (gaps between other agents' modules are not free slots)
echo "BUILD START $(date -u) tree=$P HEAD=$(git rev-parse HEAD) pkg=$PKG"
for f in $(cd $PKG/overlay && find . -type f | sed 's|^\./||'); do cmp -s $PKG/overlay/$f $P/$f || cmp -s $PKG/variants/cyclic/$f $P/$f 2>/dev/null || cmp -s $PKG/variants/packed/$f $P/$f 2>/dev/null || { echo "overlay not applied: $f"; exit 94; }; done
i=0
for m in $(grep -v '^#' $PKG/MODULES.txt | awk '{print $1}'); do
  i=$((i+1)); waitlean
  step $(printf '%02d' $i)-$m lake build $m || { echo "BUILD FAILED at $m (log $L/$(printf '%02d' $i)-$m.log)"; exit 1; }
done
waitlean
step axioms lake env lean Work/FourierC588M/Axioms.lean
step compare-fourier lake env lean --run tools/Compare.lean Work.FourierC588M.UniformFourierChallenge \
  Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
step compare-negative lake env lean --run tools/Compare.lean Work.FourierC588MNeg.UniformFourierChallenge \
  Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
NEG=$?
echo "(negative control must be nonzero: rc=$NEG)"
python3 - "$L" "$NEG" <<'PY'
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
sys.exit(0 if all(v for _, v in ok) else 1)
PY
RC=$?
grep -h "depends on axioms" $L/*C588MControls*.log $L/axioms.log 2>/dev/null
echo "BUILD END $(date -u) rc=$RC"
exit $RC
