#!/bin/bash
# i-skel: the hypothesis-mode chain L3945 (tools/mk_l3945.py hyp): install, lake build, axioms (chain + extra), greps,
# challenges re-derived from OpenAI's statement.  Compare.lean is not run (the final theorems take the engine fact as an
# instance argument by design).  usage: isk-hyp.sh <SUF> <NUM> <generated dir with stage/ and meta/>
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
SUF=$1; NUM=$2; G=$3
echo "HYP RUN $SUF NUM=$NUM START $(date -u) on $(hostname)"
[ "$(python3 -I -c "import json;print(json.load(open('$G/meta/info.json'))['mode'])")" = hyp ] || fail "not a hypothesis-mode stage"
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "stage hashes"
bash $S/isk-tree.sh || fail tree
bash $S/isk-install.sh $SUF $G/stage || fail install
bash $S/isk-build.sh $SUF $(cat $G/meta/targets.txt) || fail build
ISK_NOCOMPARE=1 bash $S/isk-checks.sh $SUF $NUM $SUF $G/meta/AxExtra$SUF.lean || fail checks
echo "HYP RUN $SUF: AS EXPECTED $(date -u)"
echo EXIT 0
