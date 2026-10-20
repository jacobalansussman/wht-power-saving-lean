#!/bin/bash
# bB-data: i-skel's checks on the BUILT final L3945H chain in my tree (axioms of Axioms.lean + AxExtra, sorry /
# forbidden greps over ENG + my 176 files + the L3945H stage, both challenges re-derived from OpenAI's statement,
# Compare.lean positive PASS / negative FAIL).  The official comparator and the bigospeedrun re-creation are the
# integration's isk-final.sh.
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
G=$R/in/final-L3945H
(cd $T && sha256sum -c --quiet $G/meta/stage.sha256) || fail "tree != final stage"
cat $R/in-ENG.files $R/in-B3945H.files $R/in-L3945Hhyp.files | sort -u > $R/in-L3945H-all.files
bash $S/isk-checks.sh L3945H 17581520 L3945H-all $G/meta/AxExtraL3945H.lean || fail "checks L3945H"
echo "CHECKS L3945H OK"
