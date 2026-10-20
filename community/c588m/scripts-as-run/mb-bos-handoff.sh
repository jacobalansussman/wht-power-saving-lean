#!/bin/bash
# m-build: take over the Lean slot of mb-final.sh when it ends (no release in between), then run the three offline
# bigospeedrun/verifier phases on the C588M stage (mb-bos-run.sh, back to back, one Lean process at a time).
L=/content/work/logs/m-build
until grep -q "FINAL END" $L/mb-final.out 2>/dev/null; do sleep 3; done
echo "HANDOFF from mb-final.sh $(date -u)"
NOWAIT=1 bash /content/work/mb/mb-bos-run.sh lean audit nanoda > $L/bos/run.log 2>&1
echo "BOS DONE $(date -u)"; grep -E "=== .* END|accepts|okay|uncaught|depends on axioms" $L/bos/run.log
