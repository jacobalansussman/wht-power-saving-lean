#!/bin/bash
# rwht-run.sh <tree root (cp -a of /content/work/shea)> <m-gen package dir .../c588m> <WHT src dir from gen_wht.py> [variant] [stop module]
# r-wht (13 Oct 2026).  On a private copy of Ryan Shea's built tree:
#   1. apply m-gen's C588M overlay (its own apply-overlay.sh; new files only, SHA256SUMS checked);
#   2. build the C588M modules ONE per `lake build`, in MODULES.txt order, up to and including the module that
#      proves the kernel engine theorem (default Work.FourierC588M.Seam) -- the WHT corollary needs nothing after it;
#   3. lay the generated WHT files over the tree (new files only) and build Challenge, Neg Challenge, Main, Solution;
#   4. checks: #print axioms, sorry/escape grep, Compare.lean (ours = published tools/Compare.lean) positive and
#      negative (one unit stronger), then a summary that prints ALL CHECKS PASSED or the failing check.
# Every Lean step runs under /usr/bin/time -v and waits while 2 or more Lean processes run on the machine.
# Swappable: everything certificate-specific comes from <src>/PARAMS.json (written by gen_wht.py) and the arguments.
set -u
ROOT=$1; PKG=$2; SRC=$3; VAR=${4:-annealed}; STOP=${5:-Work.FourierC588M.Seam}
T=$ROOT/lean/complex588
LBL=${LBL:-r-wht}; L=/content/work/logs/$LBL; mkdir -p $L
. /content/work/env.sh
export LEAN_NUM_THREADS=2
P() { python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])" $SRC/PARAMS.json $1; }
D=$(P dir); DN=$(P negdir); THM=$(P theorem); INNER=$(P inner); NUM=$(P num)
waitslot() { local w=0; while [ "$(pgrep -c -x lean)" -ge 2 ]; do sleep 10; w=$((w+10)); done; echo "  (slot free after ${w}s; lean procs now: $(pgrep -c -x lean))"; }
step() { # step <name> <command...>
  local n=$1; shift; local s=$(date +%s)
  echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log | cut -c1-200)"
  return $rc
}
T0=$(date +%s)
cd $T || exit 90
echo "RWHT START $(date -u) tree=$T HEAD=$(git rev-parse HEAD) tracked changes: $(git status --short --untracked-files=no | wc -l)"
echo "params: dir=$D negdir=$DN theorem=$THM inner=$INNER num=$NUM variant=$VAR stop=$STOP"
echo "Compare.lean: $(sha256sum /content/work/ours/Compare.lean)  tree tools/Compare.lean: $(sha256sum < tools/Compare.lean | cut -c1-64)"
echo "Work/Block/WHT.lean: $(sha256sum < Work/Block/WHT.lean | cut -c1-64)"

echo "--- 1. overlay (m-gen package $PKG, variant $VAR)"
if [ -e $L/applied.txt ] && grep -q "^variant $VAR" $L/applied.txt; then echo "already applied (applied.txt)"; else
  LBL=$LBL bash $PKG/scripts/apply-overlay.sh $T $VAR || { echo "OVERLAY FAILED"; exit 91; }; fi

echo "--- 2. C588M modules up to $STOP, one per lake call"
i=0
for m in $(grep -v '^#' $PKG/MODULES.txt | awk '{print $1}'); do
  i=$((i+1)); n=$(printf '%02d' $i)-$m
  if grep -q "=== STEP $n END rc=0" $L/run.out 2>/dev/null; then echo "  (done before: $n)"; else
    waitslot; step $n lake build $m || { echo "BUILD FAILED at $m (log $L/$n.log)"; exit 1; }; fi
  [ "$m" = "$STOP" ] && break
done

echo "--- 3. WHT overlay from $SRC (new files only)"
(cd $SRC && sha256sum -c --quiet MANIFEST.sha256) || { echo "WHT src MANIFEST mismatch"; exit 92; }
while read -r h f; do
  if [ -e "$f" ] && ! cmp -s "$SRC/$f" "$f"; then echo "COLLISION: $f exists with other bytes"; exit 93; fi
done < $SRC/MANIFEST.sha256
while read -r h f; do mkdir -p "$(dirname "$f")"; cp "$SRC/$f" "$f"; done < $SRC/MANIFEST.sha256
sha256sum -c $SRC/MANIFEST.sha256 | sed 's/^/  /'
for m in Work.$D.Challenge Work.$DN.Challenge Work.$D.Main Work.$D.Solution; do
  waitslot; step wht-build-$m lake build $m || { echo "BUILD FAILED at $m (log $L/wht-build-$m.log)"; exit 1; }
done
step wht-nobuild lake build --no-build Work.$D.Main Work.$D.Solution Work.$D.Challenge Work.$DN.Challenge $STOP

echo "--- 4. checks"
waitslot; step axioms lake env lean Work/$D/Axioms.lean; cat $L/axioms.log
echo "=== SORRY/ESCAPE GREP in Work/$D Work/$DN (Challenge files carry the statement's deliberate 'sorry')"
grep -rnwE 'sorry|native_decide|implemented_by|extern|unsafe|ofReduceBool|skipKernelTC|axiom|opaque' Work/$D Work/$DN | tee $L/sorry-grep.txt
waitslot; step compare-pos lake env lean --run /content/work/ours/Compare.lean Work.$D.Challenge Work.$D.Solution $THM; POS=$?
grep -v "  identical:" $L/compare-pos.log
waitslot; step compare-neg lake env lean --run /content/work/ours/Compare.lean Work.$DN.Challenge Work.$D.Solution $THM; NEG=$?
grep -v "  identical:" $L/compare-neg.log
echo "tracked changes in the tree after the run: $(git status --short --untracked-files=no | wc -l)"
python3 -I - "$L" "$POS" "$NEG" "$THM" "$INNER" "$D" <<'PY'
import sys, re
from pathlib import Path
L, pos, neg, thm, inner, D = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]
ax = (L/'axioms.log').read_text(); cp = (L/'compare-pos.log').read_text(); cn = (L/'compare-neg.log').read_text()
sg = (L/'sorry-grep.txt').read_text().splitlines()
std = '[propext, Classical.choice, Quot.sound]'
ok = [
 ('axioms: inner theorem on the 3 standard axioms', f"'{inner}' depends on axioms: {std}" in ax),
 ('axioms: comparator theorem on the 3 standard axioms', f"'{thm}' depends on axioms: {std}" in ax),
 ('axioms: no sorryAx', 'sorryAx' not in ax),
 ('grep: no escape hatch in Main/Solution (only comment mentions)', all(('/Challenge.lean' in l) or ('No `sorry`' in l) for l in sg)),
 ('grep: the only code sorry is the statement proof in the two Challenge files', sum(1 for l in sg if re.match(r'^[^:]+:\d+:\s+sorry\s*$', l)) == 2),
 ('Compare positive: rc 0 and RESULT: PASS', pos == 0 and 'RESULT: PASS' in cp and 'OK: statement is the same kernel term in both' in cp),
 ('Compare negative: rc != 0, RESULT: FAIL, statement differs', neg != 0 and 'RESULT: FAIL' in cn and 'FAIL: statement (type) differs' in cn),
]
for n, v in ok: print(('PASS ' if v else 'FAIL ') + n)
print('ALL CHECKS PASSED' if all(v for _, v in ok) else 'SOME CHECK FAILED')
sys.exit(0 if all(v for _, v in ok) else 1)
PY
RC=$?
echo "RWHT END $(date -u) rc=$RC total_wall=$(( $(date +%s) - T0 ))s"
exit $RC
