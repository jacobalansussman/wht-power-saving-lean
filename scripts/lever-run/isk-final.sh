#!/bin/bash
# i-skel: THE FIXED INTEGRATION RUN of a generated Fourier instance (mode final or stub of tools/mk_l3945.py) in a fresh
# private copy of ours/T: install (refuse different existing files) -> lake build of meta/targets.txt -> axioms (chain
# + meta/AxExtra<SUF>.lean) and greps -> challenges re-derived from OpenAI's statement -> Compare.lean positive PASS /
# negative FAIL -> official comparator accepted / rejected -> bigospeedrun re-creation (lean, audit, nanoda) pass ->
# bigospeedrun at one grid step too high rejected.  i-skel ran it on the stub dry run L3945S; i-final runs it on L3945.
# usage: isk-final.sh <SUF> <NUM> <generated dir with stage/ and meta/> [<extra overlay dir installed first, e.g. the
#        engine stack: Work/Mixed/..., Work/L3945/Engine.lean>] [bos threads, default 3]
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
SUF=$1; NUM=$2; G=$3; EXTRA=${4:-}; THR=${5:-3}
echo "FINAL RUN $SUF NUM=$NUM START $(date -u) on $(hostname); R=$R"
[ "$(python3 -I -c "import json;print(json.load(open('$G/meta/info.json'))['dft_numerator'])")" = "$NUM" ] || fail "meta/info.json NUM differs"
[ "$(python3 -I -c "import json;print(json.load(open('$G/meta/info.json'))['mode'])")" != hyp ] || fail "hypothesis-mode stage given"
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "stage hashes"
bash $S/isk-tree.sh || fail tree
if [ -n "$EXTRA" ] && [ "$EXTRA" != - ]; then bash $S/isk-install.sh $SUF-extra $EXTRA || fail "install extra"; fi
bash $S/isk-install.sh $SUF $G/stage || fail install
# new top-level Lean libraries of the engine stack (anything outside Work/, CSX/, XExit/): ISK_LEANLIBS="Lib1 Lib2"
for lib in ${ISK_LEANLIBS:-}; do
  grep -q "^lean_lib $lib\b" $T/lakefile.lean && echo "lakefile: lean_lib $lib present" || { printf '\nlean_lib %s where globs := #[`%s.+]\n' $lib $lib >> $T/lakefile.lean; echo "lakefile: added lean_lib $lib (globs $lib.+)"; }
done
bash $S/isk-build.sh $SUF $(cat $G/meta/targets.txt) || fail build
# greps over the stage AND the engine stack (the axioms check covers sorryAx / ofReduceBool semantically; the greps are extra)
cat $R/in-$SUF.files $( [ -e $R/in-$SUF-extra.files ] && echo $R/in-$SUF-extra.files ) | sort -u > $R/in-$SUF-all.files
bash $S/isk-checks.sh $SUF $NUM $SUF-all $G/meta/AxExtra$SUF.lean || fail checks
# the comparator (my tree's oleans) and the bigospeedrun re-creation (its own stage, compiled from the packaged sources)
# are independent: run them side by side (each waits for a heavy slot; bos uses THR Lean threads)
echo "comparator and bigospeedrun side by side: logs $LOGROOT/$SUF-cmp.out, $LOGROOT/$SUF-bos.out"
bash $S/isk-cmp.sh $SUF > $LOGROOT/$SUF-cmp.out 2>&1 &
P1=$!
( bash $S/isk-bos.sh $SUF $SUF $NUM $THR pos ${ISK_BOS_COMMENT:-} && bash $S/isk-bos.sh ${SUF}neg $SUF $NUM $THR neg ${ISK_BOS_COMMENT:-} ) > $LOGROOT/$SUF-bos.out 2>&1 &
P2=$!
wait $P1; r1=$?; wait $P2; r2=$?
grep -h VERDICT $LOGROOT/$SUF-cmp.out $LOGROOT/$SUF-bos.out
[ $r1 = 0 ] || fail "comparator (see $LOGROOT/$SUF-cmp.out)"
[ $r2 = 0 ] || fail "bigospeedrun (see $LOGROOT/$SUF-bos.out)"
ls -la $R/*-OK
echo "FINAL RUN $SUF: EVERY VERDICT AS EXPECTED $(date -u)"
echo EXIT 0
