#!/bin/bash
# m-build: after mb-build.sh passed. The official comparator (c-base's procedure: cb-cmp.sh = mb-cmp.sh, jail2.sh
# stand-in for landrun, tools at the verifier pins in /content/work/cc/bin) on the built C588M tree, positive config
# and the negative config one unit stronger; then (optional) the offline bigospeedrun/verifier phases.
# usage: mb-post.sh [cmp|bos|all]
T=/content/work/mb/shea/lean/complex588; L=/content/work/logs/m-build; WHAT=${1:-all}
waitslot() { local ok=0; while [ $ok -lt 3 ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 15; else ok=$((ok+1)); sleep 4; fi; done; }
echo "POST START $(date -u) what=$WHAT"
if [ $WHAT = cmp ] || [ $WHAT = all ]; then
  chown -R 65534:65534 $T/.lake
  waitslot; bash /content/work/mb/mb-cmp.sh m-build-A1 $T comparator/fourier-c588m.json > $L/cmp-A1.out 2>&1; echo "A1 positive rc=$? $(date -u)"
  tail -4 $L/cmp-A1.out
  waitslot; bash /content/work/mb/mb-cmp.sh m-build-A2neg $T comparator/fourier-c588m-neg.json > $L/cmp-A2neg.out 2>&1; echo "A2 negative rc=$? (must be nonzero) $(date -u)"
  grep -E "uncaught exception|does not match" $L/cmp-A2neg.out | head -3
fi
if [ $WHAT = bos ] || [ $WHAT = all ]; then
  grep -q "PREP DONE" $L/bos/prep.log || { echo "bos prep missing"; exit 3; }
  bash /content/work/mb/mb-bos-run.sh > $L/bos/run.log 2>&1; echo "bos phases done $(date -u)"; grep -E "=== .* END" $L/bos/run.log
fi
echo "POST END $(date -u)"
