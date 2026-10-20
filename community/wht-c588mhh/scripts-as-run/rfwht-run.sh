#!/bin/bash
# rfwht-run.sh (rf-wht, 13 Oct 2026, Colab #2): Walsh-Hadamard corollary at the C588MHH kernel figure.
#   0. own tree /content/work/rf-wht/shea: a copy of r-wht's tree (Shea's tree at 5239484 + the c588m overlay, annealed,
#      C588M modules built up to Work.FourierC588M.Seam), WITHOUT .lake/packages (symlink to r-wht's, read only);
#   1. overlay, new files only: the 17 h-L3-build files the corollary needs (16 of Work/XHost + Work/FourierC588MHH/Seam.lean,
#      byte-identical to h-L3-build/lean) and the eight WHT files written by gen_wht.py (MANIFEST.sha256 checked;
#      a file that exists with other bytes = stop);
#   2. compile with plain `lean` (no Lake: never run lake as root in a tree whose .lake is owned by 65534), one process
#      per module, two lanes at most (the 2-heavy-processes rule), each start waits for >= 14 GB MemAvailable;
#   3. checks: #print axioms, sorry/escape grep, Compare.lean (= published tools/Compare.lean) positive and negative.
# usage: rfwht-run.sh <stage dir with ovl/ and src-C588MHH/>      logs: /content/work/logs/rf-wht/
set -u
S=$1
R=/content/work/rf-wht/shea; T=$R/lean/complex588
SRCT=/content/work/r-wht/shea
L=/content/work/logs/rf-wht; mkdir -p $L
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
SRC=$S/src-C588MHH; OVL=$S/ovl
P() { python3 -I -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])" $SRC/PARAMS.json $1; }
D=$(P dir); DN=$(P negdir); THM=$(P theorem); INNER=$(P inner); NUM=$(P num)
T0=$(date +%s)
echo "RFWHT START $(date -u) params: dir=$D negdir=$DN theorem=$THM inner=$INNER num=$NUM"

echo "--- 0. tree copy"
if [ -e $R/.rfwht-copied ]; then echo "already copied"; else
  [ -e $R ] && { echo "$R exists but is not marked copied; stop"; exit 80; }
  mkdir -p /content/work/rf-wht
  s=$(date +%s)
  (cd /content/work/r-wht && tar --exclude=shea/lean/complex588/.lake/packages -cf - shea) | (cd /content/work/rf-wht && tar -xpf -) || { echo "COPY FAILED"; exit 81; }
  ln -s $SRCT/lean/complex588/.lake/packages $T/.lake/packages
  touch $R/.rfwht-copied
  echo "copied in $(( $(date +%s)-s ))s; .lake/packages -> $(readlink $T/.lake/packages)"
fi
cd $T || exit 82
echo "tree=$T HEAD=$(git rev-parse HEAD) tracked changes: $(git status --short --untracked-files=no | wc -l)"
echo "Compare.lean: $(sha256sum /content/work/ours/Compare.lean)  tree tools/Compare.lean: $(sha256sum < tools/Compare.lean | cut -c1-64)"
echo "Work/Block/WHT.lean: $(sha256sum < Work/Block/WHT.lean | cut -c1-64)"
echo "oleans in the copied build dir: $(find .lake/build/lib/lean -name '*.olean' | wc -l)"

echo "--- 1. overlay (new files only)"
for M in $OVL/MANIFEST.sha256 $SRC/MANIFEST.sha256; do
  B=$(dirname $M)
  (cd $B && sha256sum -c --quiet MANIFEST.sha256) || { echo "MANIFEST mismatch in $B"; exit 92; }
  while read -r h f; do
    if [ -e "$f" ] && ! cmp -s "$B/$f" "$f"; then echo "COLLISION: $f exists with other bytes"; exit 93; fi
  done < $M
  while read -r h f; do mkdir -p "$(dirname "$f")"; cp "$B/$f" "$f"; done < $M
  echo "  $(sha256sum -c $M 2>&1 | grep -c ': OK$') of $(wc -l < $M) files of $M in place"
done
echo "tracked changes after the overlay: $(git status --short --untracked-files=no | wc -l)"

LP=$T/.lake/build/lib/lean
for p in mathlib batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible Cli; do LP=$LP:$T/.lake/packages/$p/.lake/build/lib/lean; done
export LEAN_PATH=$LP LEAN_NUM_THREADS=2
# Compare.lean runs `lean` by name as a subprocess (attempt 2 failed there: "could not execute external process 'lean'")
export PATH=$(dirname $LEAN):$PATH
myheavy() { ps -C lean -o rss=,args= | awk -v t="$T/" '$1>2000000 && index($0,t)>0' | wc -l; }
memok() { [ $(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo) -ge 14 ]; }
wait_ok() { local w=0; while [ $(myheavy) -ge 2 ] || ! memok; do sleep 5; w=$((w+5)); done; [ $w -gt 0 ] && echo "  (waited ${w}s)"; return 0; }
mk() { # mk <lane> <module>...   one plain-lean compile per module; stop the lane at the first failure
  local lane=$1; shift
  for M in "$@"; do
    local f=$T/$(echo $M | tr . /).lean o=$T/.lake/build/lib/lean/$(echo $M | tr . /)
    if grep -q "^$M EXIT 0 " $L/build.summary 2>/dev/null; then echo "  (done before: $M)"; continue; fi
    wait_ok; mkdir -p $(dirname $o)
    echo "=== $lane START $M $(date -u +%H:%M:%S) MemAvailable_GB=$(awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo)"
    /usr/bin/time -v -o $L/b-$M.time $LEAN -DautoImplicit=false -R $T -o $o.olean -i $o.ilean $f > $L/b-$M.log 2>&1
    local ec=$?
    echo "$M EXIT $ec wall=$(awk -F': ' '/Elapsed/ {print $2}' $L/b-$M.time) maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/b-$M.time) errors=$(grep -c 'error' $L/b-$M.log) sorry_mentions=$(grep -ci sorry $L/b-$M.log)" | tee -a $L/build.summary
    [ $ec -ne 0 ] && { echo "BUILD FAILED at $M (log $L/b-$M.log)"; return 1; }
  done
  return 0   # (attempt 1 lacked this line: a lane whose modules all passed returned the status of the last `[ $ec -ne 0 ]`)
}

echo "--- 2. build (plain lean, LEAN_PATH = own build dir + Mathlib packages)"
mk A Work.XHost.HostK Work.XHost.HRoles Work.XHost.HDef Work.XHost.HStage Work.XHost.HCert > $L/laneA.out 2>&1 &
PA=$!
sleep 20
mk B Work.XHost.Frames Work.XHost.HSched Work.XHost.HSchedK Work.XHost.HPriceK Work.XHost.Data.H588MSched \
     Work.XHost.Data.RateC588MHHx Work.$D.Challenge Work.$DN.Challenge > $L/laneB.out 2>&1 &
PB=$!
wait $PA; RA=$?; wait $PB; RB=$?
cat $L/laneA.out $L/laneB.out
[ $RA -eq 0 ] && [ $RB -eq 0 ] || { echo "LANE FAILED (A rc=$RA, B rc=$RB)"; exit 1; }
mk C Work.XHost.HGeom Work.XHost.HEnd Work.XHost.HBank Work.XHost.Data.InstC588MHH Work.XHost.C588MHH \
     Work.FourierC588MHH.Seam Work.$D.Main Work.$D.Solution || exit 1

echo "--- 3. checks"
s=$(date +%s); wait_ok
/usr/bin/time -v -o $L/axioms.time $LEAN -DautoImplicit=false Work/$D/Axioms.lean > $L/axioms.log 2>&1; echo "axioms rc=$? wall=$(( $(date +%s)-s ))s"; cat $L/axioms.log
cat > /content/work/rf-wht/AxSeam.lean <<EOF
import Work.FourierC588MHH.Seam
#print axioms OAI.PowerSaving.SeamC588MHH.engine_C588MHHx
#print axioms OAI.PowerSaving.XHost.engine_of_ginv_host
#print axioms OAI.PowerSaving.XHost.HostSched.ofTab_fitsR
#check @OAI.PowerSaving.SeamC588MHH.engine_C588MHHx
EOF
s=$(date +%s); wait_ok
/usr/bin/time -v -o $L/axioms-seam.time $LEAN -DautoImplicit=false /content/work/rf-wht/AxSeam.lean > $L/axioms-seam.log 2>&1; echo "axioms-seam rc=$? wall=$(( $(date +%s)-s ))s"; cat $L/axioms-seam.log
echo "=== SORRY/ESCAPE GREP in Work/$D Work/$DN (Challenge files carry the statement's deliberate 'sorry')"
grep -rnwE 'sorry|native_decide|implemented_by|extern|unsafe|ofReduceBool|skipKernelTC|axiom|opaque' Work/$D Work/$DN | tee $L/sorry-grep.txt
echo "=== SORRY/ESCAPE GREP in the hosting overlay Work/XHost (without Test/) and Work/FourierC588MHH/Seam.lean"
grep -rnwE 'sorry|native_decide|implemented_by|extern|unsafe|ofReduceBool|skipKernelTC|axiom|opaque' --exclude-dir=Test Work/XHost Work/FourierC588MHH/Seam.lean | tee $L/sorry-grep-overlay.txt
cmpstep() { # cmpstep <name> <challenge module>
  local n=$1 s=$(date +%s); wait_ok
  echo "=== STEP $n START $(date -u) cmd: lean --run Compare.lean $2 Work.$D.Solution $THM"
  /usr/bin/time -v -o $L/$n.time $LEAN --run /content/work/ours/Compare.lean $2 Work.$D.Solution $THM > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time)"
  grep -v "  identical:" $L/$n.log
  return $rc
}
cmpstep compare-pos Work.$D.Challenge; POS=$?
cmpstep compare-neg Work.$DN.Challenge; NEG=$?
echo "tracked changes in the tree after the run: $(git status --short --untracked-files=no | wc -l)"
python3 -I - "$L" "$POS" "$NEG" "$THM" "$INNER" <<'PY'
import sys, re
from pathlib import Path
L, pos, neg, thm, inner = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
ax = (L/'axioms.log').read_text(); axs = (L/'axioms-seam.log').read_text()
cp = (L/'compare-pos.log').read_text(); cn = (L/'compare-neg.log').read_text()
sg = (L/'sorry-grep.txt').read_text().splitlines(); so = (L/'sorry-grep-overlay.txt').read_text().splitlines()
std = '[propext, Classical.choice, Quot.sound]'
ok = [
 ('axioms: inner theorem on the 3 standard axioms', f"'{inner}' depends on axioms: {std}" in ax),
 ('axioms: comparator theorem on the 3 standard axioms', f"'{thm}' depends on axioms: {std}" in ax),
 ('axioms: engine, hosted end theorem and checker soundness on standard axioms only, no sorryAx',
  axs.count('depends on axioms: [propext') == 3 and 'sorryAx' not in axs and 'ofReduceBool' not in axs and 'error' not in axs),
 ('axioms: no sorryAx', 'sorryAx' not in ax),
 ('grep: no escape hatch in Main/Solution (only comment mentions)', all(('/Challenge.lean' in l) or ('No `sorry`' in l) for l in sg)),
 ('grep: the only code sorry is the statement proof in the two Challenge files', sum(1 for l in sg if re.match(r'^[^:]+:\d+:\s+sorry\s*$', l)) == 2),
 ('grep: no code escape hatch in the hosting overlay (comment mentions only)', all(re.search(r'No `sorry`|no `sorry`|No sorry|`sorry`', l) for l in so)),
 ('Compare positive: rc 0 and RESULT: PASS', pos == 0 and 'RESULT: PASS' in cp and 'OK: statement is the same kernel term in both' in cp),
 ('Compare negative: rc != 0, RESULT: FAIL, statement differs', neg != 0 and 'RESULT: FAIL' in cn and 'FAIL: statement (type) differs' in cn),
]
for n, v in ok: print(('PASS ' if v else 'FAIL ') + n)
print('ALL CHECKS PASSED' if all(v for _, v in ok) else 'SOME CHECK FAILED')
sys.exit(0 if all(v for _, v in ok) else 1)
PY
RC=$?
echo "RFWHT END $(date -u) rc=$RC total_wall=$(( $(date +%s) - T0 ))s"
exit $RC
