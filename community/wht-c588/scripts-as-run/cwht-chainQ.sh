#!/bin/bash
# c-wht: after the copy of c-pack's tree finishes: build the Q modules, then axioms/grep/Compare, then the comparator.
set -u
D=/content/work/c-wht-q; L=/content/work/logs/c-wht/packed
rm -rf $D; echo "copy start $(date -u)"; s=$(date +%s)
cp -a /content/work/cpack $D || { echo "copy failed"; exit 91; }
echo "copy done $(( $(date +%s)-s ))s; chown .lake to root (a root-run Lake in a nobody-owned .lake hits git dubious-ownership and re-materialises every package)"
chown -R root:root $D/lean/complex588/.lake; echo "chown rc=$?"
bash /content/work/cwht-inQ/cwht-buildQ.sh /content/work/c-wht-q /content/work/cwht-inQ > /content/work/logs/c-wht/packed/build.log 2>&1 || { echo "build failed"; exit 92; }
CWHT_PFX=packed bash /content/work/cwht-in/cwht-checks.sh /content/work/c-wht-q OAI.PowerSaving.WHT.wht_main_block_C588Qx Work.WHTC588Q.Challenge Work.WHTC588QNeg.Challenge Work.WHTC588Q.Solution Work/WHTC588Q/Axioms.lean Work/WHTC588Q Work/WHTC588QNeg > /content/work/logs/c-wht/packed/checks.log 2>&1
# comparator: wait for the base-tree comparator runs of this agent to finish first
until grep -q '^EXIT' /content/work/logs/c-wht/cmp.log 2>/dev/null; do sleep 10; done
bash /content/work/cwht-in/cwht-cmp.sh /content/work/c-wht-q cwhtQ comparator/WHTC588Q.json comparator/WHTC588Q-nanoda.json comparator/WHTC588QNeg.json > /content/work/logs/c-wht/packed/cmp.log 2>&1
echo "chain done"
