#!/bin/bash
# Quick packed-tree re-check on the tree's own (lake-built) oleans: axioms + Compare.lean positive. Read-only.
. /content/work/env.sh
T=/content/work/c-wht-q/lean/complex588; LEAN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean
LP=""; for p in Cli batteries Qq aesop proofwidgets importGraph LeanSearchClient plausible mathlib; do LP=$LP$T/.lake/packages/$p/.lake/build/lib/lean:; done
export LEAN_PATH=$LP$T/.lake/build/lib/lean
cd $T; echo "=== START $(date -u)"; sha256sum Work/WHTC588Q/*.lean
printf "import Work.WHTC588Q.Solution\n#print axioms OAI.PowerSaving.WHTC588Q.wht_main_C588Qx\n#print axioms OAI.PowerSaving.WHT.wht_main_block_C588Qx\n#check @OAI.PowerSaving.WHT.wht_main_block_C588Qx\n" > /content/work/v-c-wht/packed-Ax.lean
/usr/bin/time -v -o /content/work/logs/v-c-wht/packed-quick-axioms.time timeout 300 nice -n 10 $LEAN /content/work/v-c-wht/packed-Ax.lean; echo "axioms rc=$?"
s=$(date +%s)
/usr/bin/time -v -o /content/work/logs/v-c-wht/packed-quick-compare.time timeout 600 nice -n 10 $LEAN --run /content/work/ours/Compare.lean Work.WHTC588Q.Challenge Work.WHTC588Q.Solution OAI.PowerSaving.WHT.wht_main_block_C588Qx 2>&1 | grep -v "  identical:"
echo "compare-pos rc=${PIPESTATUS[0]} wall=$(( $(date +%s)-s ))s"
echo "=== END $(date -u)"
