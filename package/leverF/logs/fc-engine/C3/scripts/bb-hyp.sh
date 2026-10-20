#!/bin/bash
# bB-data: the L3945H statement chain in HYPOTHESIS mode (i-skel's generator, N = 17581520) + the adapter
# Work.L3945H.EngineOfGeomW (needs Work.L3945H.Stmt and HFTest3945HW), axioms of both.
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
G=$R/in/hyp-L3945H
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "hyp stage hashes"
bash $S/isk-install.sh L3945Hhyp $G/stage || fail "install hyp"
bash $S/isk-build.sh L3945Hhyp $(cat $G/meta/targets.txt) Work.L3945H.EngineOfGeomW || fail "build hyp"
L=$LOGROOT/l3945hhyp; cd $T
wl; step axhyp lake env lean $G/meta/AxExtraL3945H.lean; cat $L/axhyp.log
axioms_ok $L/axhyp.log && ! grep -q "error" $L/axhyp.log && echo "AXIOMS hyp L3945H: standard only" || fail "axioms hyp"
grep -hE "depends on axioms" $L/build.log | grep -E "engine_of_geomW|EngineStmt" | sort -u
echo "HYP L3945H OK"
