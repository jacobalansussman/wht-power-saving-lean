#!/bin/bash
# fc-engine: controls (Ctl3945F negatives, fc-rate Rate3945FSharp) and the stub dry run L3945FS of the generator final path.
set -uo pipefail
. /content/work/fc-engine/scripts/fe-env.sh
. $ISK_SCRIPTS/isk-lib.sh
bash $S/isk-install.sh CTL $R/in/ctl || fail "install ctl"
bash $S/isk-build.sh CTL Work.Mixed.Bridge.Ctl3945F Work.Mixed.Rate3945FSharp || fail "build ctl"
grep -E "depends on axioms|error" $LOGROOT/ctl/build.log | grep -E "Ctl3945F|Rate3945FSharp|error"
G=$R/in/stub-L3945FS
(cd $G/stage && sha256sum -c --quiet ../meta/stage.sha256) || fail "stub stage hashes"
bash $S/isk-install.sh L3945FS $G/stage || fail "install stub"
bash $S/isk-build.sh L3945FS $(cat $G/meta/targets.txt) || fail "build stub"
bash $S/isk-checks.sh L3945FS 16256629 L3945FS $G/meta/AxExtraL3945FS.lean || fail "checks stub"
echo "STUB L3945FS OK"
