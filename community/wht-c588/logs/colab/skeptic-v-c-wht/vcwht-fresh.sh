#!/bin/bash
# v-c-wht: independent re-check of c-wht's claim. Compiles the overlay modules FRESH with plain `lean -o`
# into a separate directory (no lake, nothing written into the tree), then runs #print axioms and
# Compare.lean (pos + neg) against those fresh oleans, and compares fresh oleans with the tree's.
# usage: vcwht-fresh.sh <tree> <label> <Dir> <NegDir> <theorem> <inner theorem>
set -u
. /content/work/env.sh
T=$1; L=$2; D=$3; N=$4; THM=$5; INNER=$6
O=/content/work/v-c-wht/$L/olean
LOG=/content/work/logs/v-c-wht/$L
LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
mkdir -p $O/Work/$D $O/Work/$N $LOG
# Lean resolves a root package (`Work`) in the FIRST LEAN_PATH entry that has it, so $O/Work must also
# hold every other Work module: symlink each entry of the tree's build dir except the two overlay dirs.
for e in $T/.lake/build/lib/lean/Work/*; do
  b=$(basename $e); case $b in $D|$N) continue;; esac; [ -e $O/Work/$b ] || ln -s $e $O/Work/$b; done
[ -e $T/.lake/build/lib/lean/Work.olean ] && ln -sf $T/.lake/build/lib/lean/Work.olean $O/Work.olean
echo "symlinked $(find $O/Work -maxdepth 1 -type l | wc -l) entries; real dirs: $(find $O/Work -maxdepth 1 -mindepth 1 -type d | tr '\n' ' ')"
LP=$O
for p in Cli batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible mathlib; do LP=$LP:$T/.lake/packages/$p/.lake/build/lib/lean; done
LP=$LP:$T/.lake/build/lib/lean
export LEAN_PATH=$LP
cd $T || exit 90
echo "=== START $(date -u) tree=$T label=$L"
echo "git HEAD: $(git rev-parse HEAD)  status:"; git status --short
echo "LEAN_PATH=$LEAN_PATH"
$LEAN --version
sha256sum Work/$D/*.lean Work/$N/*.lean /content/work/ours/Compare.lean
waitslot() { for i in $(seq 1 18); do n=$(pgrep -c -x lean); [ "$n" -lt 3 ] && break; sleep 5; done; echo "  (lean procs running at start: $(pgrep -c -x lean))"; }
for m in $D/Challenge $N/Challenge $D/Main $D/Solution; do
  if [ -s $O/Work/$m.olean ] && [ "${REUSE_CHALLENGE:-0}" = 1 ] && [ "${m#*/}" = Challenge ]; then
    echo "--- fresh Work/$m.olean kept from attempt 1 (same source sha256, compiled $(date -u -r $O/Work/$m.olean)), not recompiled"; continue; fi
  waitslot
  s=$(date +%s)
  /usr/bin/time -v -o $LOG/fresh-$(echo $m | tr / -).time nice -n 10 $LEAN -DautoImplicit=false -R . -o $O/Work/$m.olean Work/$m.lean > $LOG/fresh-$(echo $m | tr / -).out 2>&1
  rc=$?
  echo "--- fresh lean -o Work/$m.lean rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(grep 'Maximum resident' $LOG/fresh-$(echo $m | tr / -).time | awk '{print $NF}')"
  grep -vE 'simp argument is unused|Hint:|\[apply\]|linter|^\s' $LOG/fresh-$(echo $m | tr / -).out | head -20
done
echo "--- fresh olean sha256 vs the tree's olean (built by c-wht's lake build)"
for m in $D/Challenge $N/Challenge $D/Main $D/Solution; do
  a=$(sha256sum < $O/Work/$m.olean | cut -c1-16); b=$(sha256sum < .lake/build/lib/lean/Work/$m.olean | cut -c1-16)
  echo "  $m fresh=$a tree=$b $( [ "$a" = "$b" ] && echo IDENTICAL || echo DIFFERENT)"
done
waitslot
echo "--- axioms (fresh oleans first on LEAN_PATH)"
printf 'import Work.%s.Solution\n#print axioms %s\n#print axioms %s\n#check @%s\n' $D $INNER $THM $THM > /content/work/v-c-wht/$L/Ax.lean
/usr/bin/time -v -o $LOG/axioms.time nice -n 10 $LEAN /content/work/v-c-wht/$L/Ax.lean; echo "axioms rc=$?"
for neg in pos neg; do
  waitslot
  C=Work.$D.Challenge; [ $neg = neg ] && C=Work.$N.Challenge
  echo "--- Compare.lean $neg: $C vs Work.$D.Solution $THM"
  s=$(date +%s)
  /usr/bin/time -v -o $LOG/compare-$neg.time nice -n 10 $LEAN --run /content/work/ours/Compare.lean $C Work.$D.Solution $THM 2>&1 | grep -v "  identical:"
  echo "compare-$neg rc=${PIPESTATUS[0]} wall=$(( $(date +%s)-s ))s maxrss_kb=$(grep 'Maximum resident' $LOG/compare-$neg.time | awk '{print $NF}')"
done
echo "=== END $(date -u)"
