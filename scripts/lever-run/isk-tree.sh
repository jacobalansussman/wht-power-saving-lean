#!/bin/bash
# i-skel: my private copy of our built stack (ours/T, U8X0HQ prebuilt), as COLABS.md prescribes.  Idempotent.
# usage: isk-tree.sh            (uses ISK_R / ISK_TREE, default /content/work/i-skel/T)
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
L=$LOGROOT/tree; mkdir -p $L
[ -e /content/work/OURS-READY ] || fail "no /content/work/OURS-READY on this Colab"
if [ ! -e $TT ]; then s=$(date +%s); mkdir -p $(dirname $TT); cp -a /content/work/ours/T $TT || fail "cp -a"; echo "copy of ours/T done in $(( $(date +%s)-s ))s"; fi
cd $T || fail "no $T"
if [ "$(stat -c %U .lake)" != root ]; then echo "chown .lake back to root (after a comparator run)"; chown -R root:root .lake; fi
[ "$(stat -c %U .lake)" = root ] || fail ".lake not root-owned"
step nobuild lake build --no-build Work.FourierU8X0HQ.Main || fail "ours/T copy not up to date"
tail -1 $L/nobuild.log
[ "$(sha256sum tools/Compare.lean | cut -c1-64)" = $COMPARE_SHA ] || fail "tools/Compare.lean is not e48d6d3b"
echo "TREE OK $TT ($(du -sh $TT | cut -f1))"
