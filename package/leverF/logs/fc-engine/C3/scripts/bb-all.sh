#!/bin/bash
# bB-data: the whole C3 sequence (one lake at a time in my tree): wait for the pre-build, then build, hyp, neg,
# engine (+ final chain), stub.  Each step's verdict line in /content/work/logs/bB-data/all.out.
S=/content/work/bB-data/scripts; LG=/content/work/logs/bB-data
while tmux -L bB-data has-session -t pre 2>/dev/null; do sleep 15; done
# (pre: every target built except Work.Mixed.SeamConst, which is not in ENG: it comes with the L3945H stage)
grep -q "Some required targets logged failures" $LG/pre.out && { grep -A3 "Some required targets logged failures" $LG/pre.out | grep -q "^- Work.Mixed.SeamConst$" && [ "$(grep -A9 "Some required targets logged failures" $LG/pre.out | grep -c "^- ")" = 1 ] || { echo "PRE FAILED beyond SeamConst"; exit 1; }; }
echo "pre done $(date -u +%T)"
bash $S/bb-build.sh > $LG/build.out 2>&1; tail -1 $LG/build.out; grep -q "BUILD 3945H OK" $LG/build.out || exit 1
bash $S/bb-hyp.sh > $LG/hyp.out 2>&1; tail -1 $LG/hyp.out
bash $S/bb-neg.sh > $LG/neg.out 2>&1; tail -1 $LG/neg.out
grep -q "HYP L3945H OK" $LG/hyp.out && { bash $S/bb-engine.sh > $LG/engine.out 2>&1; tail -1 $LG/engine.out; }
bash $S/bb-stub.sh > $LG/stub.out 2>&1; tail -1 $LG/stub.out
echo "ALL STEPS DONE $(date -u +%T)"
