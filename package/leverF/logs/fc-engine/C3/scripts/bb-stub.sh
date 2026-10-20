#!/bin/bash
# bB-data: stub dry run of the generator's final path with suffix L3945HS (U8X0HQ's engine as the stub, N = 16256629):
# install, build, then isk-checks.sh (axioms, greps, challenges re-derived, Compare.lean positive PASS / negative FAIL).
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
G=$R/in/stub-L3945HS
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "stub stage hashes"
bash $S/isk-install.sh L3945HS $G/stage || fail "install stub"
bash $S/isk-build.sh L3945HS $(cat $G/meta/targets.txt) || fail "build stub"
bash $S/isk-checks.sh L3945HS 16256629 L3945HS $G/meta/AxExtraL3945HS.lean || fail "checks stub"
echo "STUB L3945HS OK"
