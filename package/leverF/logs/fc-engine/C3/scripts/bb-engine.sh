#!/bin/bash
# bB-data: the integration draft Work.L3945H.Engine (Prov + EngineOfGeomW on design B's schedule, p-wgen's two words),
# then, if it compiles, the FINAL-mode L3945H chain in place of the hypothesis-mode one (the generated final stage's
# downstream files overwrite the hypothesis ones; hashes checked against the final stage).
set -uo pipefail
export ISK_R=/content/work/bB-data ISK_LOG=/content/work/logs/bB-data ISK_SCRIPTS=/content/work/bB-data/scripts
. $ISK_SCRIPTS/isk-lib.sh
bash $S/isk-build.sh L3945Hengine Work.L3945H.Engine || fail "engine"
G=$R/in/final-L3945H
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "final stage hashes"
cd $T; n=0
while read -r h f; do if ! cmp -s "$G/stage/$f" "$T/$f"; then cp -a "$G/stage/$f" "$T/$f"; n=$((n+1)); echo "final over hyp: $f"; fi; done < $G/meta/stage.sha256
echo "replaced $n files by the final stage"
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) && (cd $T && sha256sum -c --quiet $G/meta/stage.sha256) || fail "tree != final stage"
bash $S/isk-build.sh L3945Hfinal $(cat $G/meta/targets.txt) || fail "final build"
L=$LOGROOT/l3945hfinal; cd $T
wl; step axfinal lake env lean $G/meta/AxExtraL3945H.lean; cat $L/axfinal.log
axioms_ok $L/axfinal.log && ! grep -q "error" $L/axfinal.log && echo "AXIOMS final L3945H: standard only" || fail "axioms final"
echo "FINAL CHAIN L3945H BUILT"
