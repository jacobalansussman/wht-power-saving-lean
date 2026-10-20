# i-skel (round wht45/lever, 19 Oct 2026): shared helpers of the integration scripts (sourced, not run).
# Paths default to i-skel's own; i-final sets ISK_R / ISK_LOG / ISK_TREE to its own folder (see RUNBOOK.md).
# Recipes copied from fsp-build (colab/fspb-*.sh) and e-r45m1 (tools/mxt-build.sh); only paths and parameters changed.
. /content/work/env.sh
R=${ISK_R:-/content/work/i-skel}                 # my folder on the Colab
LOGROOT=${ISK_LOG:-/content/work/logs/i-skel}     # my logs
TT=${ISK_TREE:-$R/T}; T=$TT/lean/complex588       # my private copy of ours/T (Lake project)
S=${ISK_SCRIPTS:-$R/scripts}
REF=/content/work/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean
REF_SHA=$(sha256sum $REF | cut -c1-64)
TC=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin
J=/content/work/cc-jail                           # prep's jail kit (read only)
V=/content/work/bos/bos-verifier                  # clean bigospeedrun/verifier clone at c92b9b6 (read only)
COMPARE_SHA=e48d6d3b5be7f7c84ea949f1dace0a95bd58710c7cd5f3b94c60f4785975dae9
MAXHEAVY=${ISK_MAXHEAVY:-3}
mkdir -p $LOGROOT
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
avail() { free -m | awk '/Mem:/ {print $7}'; }
# wait for a heavy-Lean slot: fewer than MAXHEAVY heavy Lean processes and at least $1 MB available
wl() { local need=${1:-14000} w=0; while [ "$(heavy)" -ge $MAXHEAVY ] || [ "$(avail)" -lt $need ]; do sleep 10; w=$((w+10)); done
  echo "   (slot: waited ${w}s; heavy $(heavy), avail $(avail) MB)"; }
# step LABEL cmd...: run with /usr/bin/time into $L/LABEL.log / .time
step() { local lab=$1; shift; local s=$(date +%s); echo "=== $lab START $(date -u +%T)"
  /usr/bin/time -v -o $L/$lab.time "$@" > $L/$lab.log 2>&1; local rc=$?
  echo "=== $lab END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$lab.time) $(date -u +%T)"; return $rc; }
fail() { echo "FAILED: $*"; echo "EXIT 1"; exit 1; }
# axioms_ok LOG: every '... depends on axioms: [...]' line lists only propext / Classical.choice / Quot.sound,
# no 'does not depend' surprises are failures; at least one line must be present.
axioms_ok() { local f=$1
  grep -q "depends on axioms" $f || grep -q "does not depend on any axioms" $f || return 1
  ! grep -E "depends on axioms" $f | sed -E 's/.*\[(.*)\].*/\1/' | tr ',' '\n' | sed 's/ //g' | grep -v -x -E 'propext|Classical\.choice|Quot\.sound' | grep -q . ; }
lower() { echo "$1" | tr A-Z a-z; }
