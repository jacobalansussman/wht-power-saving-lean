#!/bin/bash
# rf-fresh2.sh <package name, e.g. leverF> : rf-fresh.sh with three changes (the bundle is a setting, BUNDLE=...;
# SKIP_U8=1 leaves out the optional package/u8x0hq step; MIN_FREE_MB makes the Compare and comparator steps wait
# for free memory). The script as it was run on the Colab machines; its paths are theirs.
set -u
P=${1:?package}; W=/content/work/r-finalize; LOG=/content/work/logs/r-finalize/$P; mkdir -p $LOG; BUNDLE=${BUNDLE:-$W/in/r7test.bundle}
. /content/work/env.sh
echo "START $(date -u) package=$P"
rm -rf $W/repo-$P $W/shea-$P
git clone -q -b main $BUNDLE $W/repo-$P || { echo "CLONE FAILED"; exit 1; }
git -C $W/repo-$P log --oneline -2; echo "tree $(git -C $W/repo-$P rev-parse HEAD^{tree})"
git -C $W/repo-$P status --short | wc -l
s=$(date +%s); cp -a /content/work/shea $W/shea-$P; echo "base copy $(( $(date +%s)-s ))s HEAD $(git -C $W/shea-$P rev-parse HEAD) status lines $(git -C $W/shea-$P status --short | wc -l)"
T=$W/shea-$P/lean/complex588
cd $W/repo-$P
echo "=== apply-overlay"; LOGDIR=$LOG/apply bash scripts/apply-overlay.sh package/$P $T; echo "apply rc=$?"
M=$(grep -v '^#' package/$P/MODULES.txt | awk 'NF {print $1}')
echo "=== lake build of all modules at once (VERIFY 12.3)"; s=$(date +%s)
( cd $T && LEAN_NUM_THREADS=2 /usr/bin/time -v -o $LOG/bulk-build.time lake build $M > $LOG/bulk-build.log 2>&1 ); echo "bulk build rc=$? $(( $(date +%s)-s ))s; $(tail -1 $LOG/bulk-build.log)"
grep -E "Maximum resident" $LOG/bulk-build.time
echo "=== build-and-check.sh package/$P"; s=$(date +%s)
LOGDIR=$LOG/bac MAX_LEAN=2 MIN_FREE_MB=${MIN_FREE_MB:-0} bash scripts/build-and-check.sh package/$P $T > $LOG/bac.out 2>&1; echo "build-and-check rc=$? $(( $(date +%s)-s ))s"
grep -E "^(PASS|FAIL|ALL|SOME|   hit)" $LOG/bac.out
if [ -z "${SKIP_U8:-}" ]; then
echo "=== build-and-check.sh package/u8x0hq (optional step of 12.3)"; s=$(date +%s)
LOGDIR=$LOG/bac-u8x0hq MAX_LEAN=2 bash scripts/build-and-check.sh package/u8x0hq $T > $LOG/bac-u8x0hq.out 2>&1; echo "build-and-check u8x0hq rc=$? $(( $(date +%s)-s ))s"
grep -E "^(ALL|SOME|FAIL)" $LOG/bac-u8x0hq.out
else echo "=== build-and-check.sh package/u8x0hq skipped (SKIP_U8)"; fi
POS=$(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d.get("comparator",{}); print(c["positive"].removeprefix("overlay/"))' package/$P/result.json)
NEG=$(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d.get("comparator",{}); print(c["negative"].removeprefix("overlay/"))' package/$P/result.json)
echo "configs from result.json: pos=$POS neg=$NEG"
memwait() { while [ "$(free -m | awk '/Mem:/ {print $7}')" -lt "${MIN_FREE_MB:-0}" ]; do sleep 30; done; }
memwait; echo "=== run-comparator.sh neg"; LOGDIR=$LOG/cmp CC=/content/work/cc JAIL=/content/work/cc-jail bash scripts/run-comparator.sh neg $T $NEG > $LOG/cmp-neg.out 2>&1; echo "comparator neg rc=$? (must be nonzero)"; grep -E "Const does not match|accepts|okay|rror" $LOG/cmp-neg.out | head -5
memwait; echo "=== run-comparator.sh pos"; LOGDIR=$LOG/cmp CC=/content/work/cc JAIL=/content/work/cc-jail bash scripts/run-comparator.sh pos $T $POS > $LOG/cmp-pos.out 2>&1; echo "comparator pos rc=$?"; grep -E "accepts|okay|rror|Const" $LOG/cmp-pos.out | head -8
grep -E "Maximum resident|Elapsed" $LOG/cmp/*.time
echo "RF-FRESH DONE $(date -u)"
