#!/bin/bash
# v-fs-statement (skeptic, 14 Oct 2026 UTC): independent read-only re-checks of fs-build's U8XMHH tree on Colab #2.
# Plain lean (no Lake; writes no olean), as uid 65534 (owner of .lake; oleans are 600); one heavy Lean process at a time.
set -uo pipefail
V=/content/work/v-fs-statement; T=/content/work/fs-build/B/lean/complex588; O=$V/out
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
mkdir -p $O; chown -R 65534:65534 $V
echo "VFS START $(date -u)"; echo "toolchain file: $(cat $T/lean-toolchain)"
LP=$T/.lake/build/lib/lean; for p in $T/.lake/packages/*; do LP=$LP:$p/.lake/build/lib/lean; done
echo "LEAN_PATH=$LP"
cd $T
echo "--- source hashes in the tree"
sha256sum Work/FourierU8XMHH/*.lean Work/FourierU8XMHHNeg/*.lean Work/XHost/Data/HU8XMSched.lean Work/XHost/Data/InstU8XMHH.lean Work/XHost/Data/RateU8XMHHx.lean Work/XHost/U8XMHH.lean tools/Compare.lean > $O/tree-hashes.txt; wc -l < $O/tree-hashes.txt
echo "--- olean vs source mtimes (olean must be newer)"
for m in FourierU8XMHH/UniformFourierChallenge FourierU8XMHHNeg/UniformFourierChallenge FourierU8XMHH/Main FourierU8XMHH/Goal FourierU8XMHH/Seam XHost/U8XMHH XHost/Data/HU8XMSched XHost/Data/RateU8XMHHx; do
  s=$(stat -c %Y Work/$m.lean); o=$(stat -c %Y .lake/build/lib/lean/Work/$m.olean); echo "$m src=$(date -u -d @$s +%T) olean=$(date -u -d @$o +%T) $([ $o -ge $s ] && echo ok || echo STALE)"; done
run() { local lab=$1; shift; while [ "$(ps -C lean -o rss= | awk '$1>2000000' | wc -l)" -ge 1 ] || [ "$(free -m | awk '/Mem:/ {print $7}')" -lt 20000 ]; do sleep 10; done
  local s=$(date +%s); echo "=== $lab START $(date -u +%T)"
  env -i PATH=$TC:/usr/bin:/bin HOME=$V LANG=C.UTF-8 LEAN_PATH=$LP LEAN_NUM_THREADS=1 /usr/bin/time -v -o $O/$lab.time setpriv --reuid 65534 --regid 65534 --clear-groups -- "$@" > $O/$lab.log 2>&1; local rc=$?
  echo "=== $lab END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $O/$lab.time)"; return 0; }
run check lean $V/VfsCheck.lean; cat $O/check.log
run cmp-pos lean --run tools/Compare.lean Work.FourierU8XMHH.UniformFourierChallenge Work.FourierU8XMHH.Main OAI.PowerSaving.transform_mainU8XMHH OAI.PowerSaving.convolution_mainU8XMHH; tail -15 $O/cmp-pos.log
run cmp-neg lean --run tools/Compare.lean Work.FourierU8XMHHNeg.UniformFourierChallenge Work.FourierU8XMHH.Main OAI.PowerSaving.transform_mainU8XMHH OAI.PowerSaving.convolution_mainU8XMHH; grep -E "FAIL|RESULT|rror" $O/cmp-neg.log | head -8
echo "--- unwritable check: any file in the tree newer than VFS start?"; find $T -newer $V/VfsCheck.lean -type f 2>/dev/null | head -5
echo "VFS END $(date -u)"
