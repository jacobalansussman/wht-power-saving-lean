#!/bin/bash
# h-L3-build: official comparator (c-base/m-build procedure, lake env, jail2.sh) positive then negative, one at a time.
N=$1; lc=$(echo $N | tr A-Z a-z); T=/content/work/hl3/shea/lean/complex588; L=/content/work/logs/h-L3-build
waitslot() { while [ "$(ps -C lean -o rss= | awk "\$1>2000000" | wc -l)" -ge 2 ]; do sleep 5; done; }
waitslot; bash /content/work/hl3/hl3-cmp.sh hl3-$N-pos $T comparator/fourier-$lc.json > $L/final-$N-cmp-pos.out 2>&1; echo "COMPARATOR positive rc=$? $(date -u +%H:%M:%S)"
grep -E "accepts the solution|Your solution is okay|uncaught exception|=== .* END" $L/final-$N-cmp-pos.out
waitslot; bash /content/work/hl3/hl3-cmp.sh hl3-$N-neg $T comparator/fourier-$lc-neg.json > $L/final-$N-cmp-neg.out 2>&1; echo "COMPARATOR negative rc=$? (must be nonzero) $(date -u +%H:%M:%S)"
grep -E "accepts the solution|Your solution is okay|uncaught exception|=== .* END|differ" $L/final-$N-cmp-neg.out | head -5
echo CMP DONE
