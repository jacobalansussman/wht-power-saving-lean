#!/bin/bash
# x-words (wht45/lever2979, 19 Oct 2026): p-gates' gates (unchanged tools) on the two #2979 words lowered by
# p-lower's tools. Light Python. Exits nonzero on any failure. W14 (price) is x-gate's; here only W0-W13, the
# second implementation, the native gate (on OUR transport of the BORROW-WORD: #2979 ships no native file) and the
# nine negatives.
set -euo pipefail
L=${LEVER:-$HOME/misc/fft/verify/checks/wht45/lever}
X=${XW:-$HOME/misc/fft/verify/checks/wht45/lever2979/x-words}
T=$L/p-gates/tools; O=$X/out/gates; PL=$X/out/lower
D=${D2979:-$L/p-alt/dl/main2979/research/rank7-joint-endframes/data}
PY="python3 -I -B"
mkdir -p $O
for x in base extra; do
  pv=$([ $x = base ] && echo gauge44 || echo extra)
  $PY $T/gw.py $PL/$x.gcert1w.json.gz $O/$x.gw.json --native $PL/$x.canonical.json.gz --word $D/WORD-$x.json.gz --label "#2979 $x ($pv, gcert/1w)"
  $PY $T/gwmirror.py $PL/$x.gcert1w.json.gz $O/$x.mirror.json --label "#2979 $x" --compare 1
  $PY $T/gwnative.py $PL/$x.canonical.json.gz $O/$x.native.json --label "#2979 $x native (our transport)" --negatives
  $PY $T/gwneg.py $PL/$x.gcert1w.json.gz $O/negatives-$x.json
done
echo GATES2979-ALL-PASS
