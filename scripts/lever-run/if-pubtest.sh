#!/bin/bash
# i-final (wht45/lever): the release candidate's PUBLISHED scripts on package/lever (INTEGRATION step 7), Colab C5
# (home), on a fresh copy of Shea's built base: apply-overlay.sh (lays package/u8x0hq first, as result.json
# "requires"), the laid tree against the hash list of the tree that was checked on C4, one parallel `lake build` of
# every module (so that the script's per-module steps replay instead of building one module at a time), then
# build-and-check.sh exactly as published.  Read only on /content/work/shea; my own copy only.
set -u
. /content/work/env.sh
W=/content/work/i-final; P=$W/pub; L=/content/work/logs/i-final/pubtest; mkdir -p $L
heavy() { ps -C lean -o rss= 2>/dev/null | awk '$1>2000000' | wc -l; }
avail() { free -m | awk '/Mem:/ {print $7}'; }
echo "PUBTEST START $(date -u) host=$(hostname) heavy=$(heavy) avail=$(avail)"
[ -e $W/shea ] && { echo "refusing: $W/shea exists"; exit 1; }
s=$(date +%s); cp -a /content/work/shea $W/shea || exit 2; echo "copy $(( $(date +%s)-s ))s $(du -sh $W/shea | cut -f1)"
T=$W/shea/lean/complex588
[ "$(stat -c %U $T/.lake)" = root ] || { echo ".lake not root"; exit 3; }
echo "base HEAD $(git -C $W/shea rev-parse HEAD); tracked changes: $(git -C $W/shea status --porcelain --untracked-files=no | wc -l)"
cd $P
LOGDIR=$L/apply bash in/scripts/apply-overlay.sh in/package/lever $T; echo "apply rc=$?"
tail -2 $L/apply/applied.txt; tail -2 $L/apply/requires-u8x0hq/applied.txt
(cd $T && find . -path ./.lake -prune -o -type f -print | sed 's#^\./##' | LC_ALL=C sort | xargs -d '\n' sha256sum) > $L/laid-tree.sha256
if diff <(cat in/SHA256SUMS.tree) $L/laid-tree.sha256 > $L/laid-vs-checked.diff; then echo "LAID TREE = CHECKED TREE: all $(wc -l < $L/laid-tree.sha256) files byte-identical (lakefile.lean included)"; else echo "LAID TREE DIFFERS from the checked tree:"; head -20 $L/laid-vs-checked.diff; fi
while [ "$(heavy)" -ge 1 ] || [ "$(avail)" -lt 40000 ]; do sleep 30; done
MODS=$(grep -v '^#' in/package/lever/MODULES.txt | awk 'NF {print $1}' | tr '\n' ' ')
U8=$(python3 -I -c 'import json;print(" ".join(json.load(open("in/package/u8x0hq/result.json"))["build_targets"]))')
s=$(date +%s); echo "=== prebuild START $(date -u)"
(cd $T && /usr/bin/time -v -o $L/prebuild.time env LEAN_NUM_THREADS=3 lake build $U8 $MODS > $L/prebuild.log 2>&1); rc=$?
echo "=== prebuild END rc=$rc wall=$(( $(date +%s)-s ))s $(grep 'Maximum resident' $L/prebuild.time | sed 's/^\s*//')"; tail -2 $L/prebuild.log
[ $rc = 0 ] || { grep -E "✖|error" $L/prebuild.log | head -20; echo "PUBTEST END (prebuild failed) $(date -u)"; exit 4; }
LOGDIR=$L/check MAX_LEAN=3 MIN_FREE_MB=20000 bash in/scripts/build-and-check.sh in/package/lever $T 2>&1 | tail -40; echo "check rc=${PIPESTATUS[0]}"
echo "PUBTEST END $(date -u)"
