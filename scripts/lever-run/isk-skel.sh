#!/bin/bash
# i-skel on C5: (1) the hypothesis-mode L3945 chain, (2) the stub dry run L3945S through the fixed final run, in ONE
# private tree (the comparator chowns .lake, so the hypothesis build goes first).
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
bash $S/isk-hyp.sh L3945 17501109 $R/in/hyp-L3945 > $LOGROOT/hyp.out 2>&1; echo "hyp: $(tail -1 $LOGROOT/hyp.out)"
bash $S/isk-final.sh L3945S 16256629 $R/in/stub-L3945S - ${1:-3} > $LOGROOT/stub.out 2>&1; echo "stub: $(tail -1 $LOGROOT/stub.out)"
echo EXIT 0
