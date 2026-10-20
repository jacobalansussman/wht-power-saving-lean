#!/bin/bash
# fc-engine: FINAL-mode L3945F chain (full cut on design B) in my copy of bB-data tree, then i-skel checks.
set -uo pipefail
. /content/work/fc-engine/scripts/fe-env.sh
. $ISK_SCRIPTS/isk-lib.sh
G=$R/in/final-L3945F
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "final stage hashes"
(cd $T && sha256sum -c --quiet $G/meta/stage.sha256) || fail "tree != final stage"
cat /content/work/bB-data/in-ENG.files /content/work/bB-data/in-B3945H.files $R/in-ENG.files $R/in-L3945F.files | sort -u > $R/in-L3945F-all.files
echo "files for the greps: $(wc -l < $R/in-L3945F-all.files)"
bash $S/isk-build.sh L3945Ffinal $(cat $G/meta/targets.txt) || fail "final build"
bash $S/isk-checks.sh L3945F 17806935 L3945F-all $G/meta/AxExtraL3945F.lean || fail "checks L3945F"
echo "FINAL L3945F BUILT AND CHECKED"
