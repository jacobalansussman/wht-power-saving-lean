#!/bin/bash
# c-wht: axioms, sorry grep, Compare.lean (positive and one-unit-stronger negative control).
# usage: cwht-checks.sh <tree root> <theorem> <challenge> <neg challenge> <solution> <axioms file> <newdirs...>
set -u
T=$1/lean/complex588; TH=$2; CH=$3; NEG=$4; SOL=$5; AX=$6; shift 6
L=/content/work/logs/c-wht/${CWHT_PFX:-base}; mkdir -p $L
. /content/work/env.sh
cd "$T" || exit 90
waitslot() { local w=0; while [ $(pgrep -c -x lean) -ge 2 ]; do sleep 10; w=$((w+10)); done; echo "  (slot free after waiting ${w}s; lean procs now: $(pgrep -c -x lean))"; }
echo "=== AXIOMS $(date -u)"
waitslot
/usr/bin/time -v -o $L/axioms.time lake env lean $AX; echo "rc=$?"; grep -E 'Maximum resident|Elapsed' $L/axioms.time
echo "=== SORRY/ESCAPE GREP in the new modules: $*"
for d in "$@"; do grep -rnwE 'sorry|native_decide|implemented_by|extern|unsafe|ofReduceBool|skipKernelTC|axiom|opaque' $d || echo "  (none in $d)"; done
echo "=== COMPARE positive $(date -u)  lean procs running: $(pgrep -c lean)"
waitslot
s=$(date +%s)
/usr/bin/time -v -o $L/compare-pos.time lake env lean --run /content/work/ours/Compare.lean $CH $SOL $TH; echo "rc=$? wall=$(( $(date +%s)-s ))s"
grep -E 'Maximum resident|Elapsed|Exit status' $L/compare-pos.time
echo "=== COMPARE negative $(date -u)  lean procs running: $(pgrep -c lean)"
waitslot
s=$(date +%s)
/usr/bin/time -v -o $L/compare-neg.time lake env lean --run /content/work/ours/Compare.lean $NEG $SOL $TH; echo "rc=$? wall=$(( $(date +%s)-s ))s"
grep -E 'Maximum resident|Elapsed|Exit status' $L/compare-neg.time
echo "=== END $(date -u)"
