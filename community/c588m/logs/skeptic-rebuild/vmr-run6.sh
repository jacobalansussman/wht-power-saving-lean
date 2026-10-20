#!/bin/bash
# v-m-rebuild (13 Oct 2026), no tree copy.  Sources = m-gen's PACKAGE overlay (checked against SHA256SUMS), compiled with
# plain `lean` into my own output roots; LEAN_PATH resolves everything else to read-only oleans.  Nothing is written into
# any other agent's tree or into c-base's base tree.
#  PHASE 1: Work.FourierC588M.Main (the final module) recompiled from package source against m-build's oleans of its
#           imports; then Deps.lean (proof-closure probe), Compare.lean positive + negative using MY Main.olean.
#  PHASE 2: all 70 overlay modules recompiled from package source against c-base's UNTOUCHED base oleans only
#           (no m-build olean on the path); then Axioms, Deps, Compare pos/neg; olean comparison with m-build's.
. /content/work/env.sh
set -u
W=/content/work/v-m-rebuild; L=/content/work/logs/v-m-rebuild; mkdir -p $L
PKG=${PKG:-/content/work/m-gen/c588m}; B=/content/work/shea/lean/complex588; MB=${MB:-/content/work/mb/shea/lean/complex588}
SRC=$W/src; O1=$W/out1; O2=$W/out2
export LEAN_NUM_THREADS=2; LEANBIN=/content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lean; [ -x $LEANBIN ] || { echo "no $LEANBIN"; exit 89; }
step() { local n=$1; shift; local s=$(date +%s); echo "=== STEP $n START $(date -u) cmd: $*"
  /usr/bin/time -v -o $L/$n.time "$@" > $L/$n.log 2>&1; local rc=$?
  echo "=== STEP $n END rc=$rc wall=$(( $(date +%s) - s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$n.time) last: $(tail -1 $L/$n.log)"; return $rc; }
waitlean() { local ok=0; while [ $ok -lt ${WAITN:-3} ]; do if [ "$(pgrep -c lean)" -ge 2 ]; then ok=0; sleep 10; else ok=$((ok+1)); sleep ${WAITS:-3}; fi; done; }
pkgpaths() { local p; for p in $B/.lake/packages/*/; do printf ':%s.lake/build/lib/lean' "$p"; done; }
MODS=$(grep -v '^#' $PKG/MODULES.txt | awk '{print $1}')
echo "START $(date -u)"
rm -rf $SRC $O1 $O2; mkdir -p $SRC $O1 $O2
(cd $PKG && sha256sum -c --quiet SHA256SUMS) || { echo "package SHA256SUMS FAIL"; exit 91; }
cp -R $PKG/overlay/Work $SRC/
mm=0; n=0; while read h f; do case "$f" in overlay/Work/*) r=${f#overlay/}; n=$((n+1)); t=$(sha256sum "$SRC/$r" | cut -d' ' -f1); [ "$t" = "$h" ] || { echo "HASH MISMATCH $r"; mm=1; }
  if [ -d $MB ]; then t2=$(sha256sum "$MB/$r" | cut -d' ' -f1); [ "$t2" = "$h" ] || { echo "M-BUILD TREE DIFFERS $r"; mm=1; }; fi;; esac; done < $PKG/SHA256SUMS
echo "sources: $n .lean files checked against SHA256SUMS (my copy and m-build's tree); mismatches: $mm"
echo "modules in MODULES.txt: $(echo $MODS | wc -w)"
pre=0; for m in $MODS; do f=$(echo $m | tr . /); [ -e $B/.lake/build/lib/lean/$f.olean ] && { echo "BASE HAS $f.olean"; pre=1; }; done
echo "overlay modules already present as oleans in c-base's base tree: $pre (must be 0)"
## PHASE 1
if [ "${PHASE1:-1}" = 1 ]; then
cp -rs $MB/.lake/build/lib/lean/. $O1/ || exit 83
LP1=$O1$(pkgpaths)
f=Work/FourierC588M/Main; [ -L $O1/$f.olean ] || exit 84; rm -f $O1/$f.olean $O1/$f.ilean $O1/$f.olean.hash $O1/$f.ilean.hash
[ -e $O1/$f.olean ] && exit 85; echo "P1 farm: $(find $O1 -type l | wc -l) symlinks to m-build oleans; Main removed from the farm"
waitlean; step P1-main env LEAN_PATH=$LP1 $LEANBIN -DautoImplicit=false -R $SRC $SRC/$f.lean -o $O1/$f.olean -i $O1/$f.ilean
cmp -s $O1/$f.olean $MB/.lake/build/lib/lean/$f.olean && echo "P1: my Main.olean is byte-identical to m-build's" || echo "P1: my Main.olean differs in bytes from m-build's"
waitlean; step P1-deps env LEAN_PATH=$LP1 $LEANBIN -DautoImplicit=false $W/Deps.lean
waitlean; step P1-compare-pos env LEAN_PATH=$LP1 $LEANBIN --run /content/work/ours/Compare.lean Work.FourierC588M.UniformFourierChallenge Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
waitlean; step P1-compare-neg env LEAN_PATH=$LP1 $LEANBIN --run /content/work/ours/Compare.lean Work.FourierC588MNeg.UniformFourierChallenge Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
echo "P1 negative rc=$? (must be nonzero)"
echo "PHASE1 END $(date -u)"
fi
[ "${PHASE2:-1}" = 1 ] || { echo "END (phase 2 skipped) $(date -u)"; exit 0; }
## PHASE 2
cp -rs $B/.lake/build/lib/lean/. $O2/ || exit 86
LP2=$O2$(pkgpaths)
echo "P2 farm: $(find $O2 -type l | wc -l) symlinks to c-base base oleans, $(find $O2 -type f | wc -l) regular files"
i=0
for m in $MODS; do
  i=$((i+1)); f=$(echo $m | tr . /); mkdir -p $O2/$(dirname $f); { [ -e $O2/$f.olean ] || [ -L $O2/$f.olean ] || [ -L $O2/$f.ilean ]; } && { echo "REFUSE: $f exists in farm"; exit 87; }
  waitlean; step P2-$(printf '%02d' $i)-$m env LEAN_PATH=$LP2 $LEANBIN -DautoImplicit=false -R $SRC $SRC/$f.lean -o $O2/$f.olean -i $O2/$f.ilean || { echo "P2 BUILD FAILED at $m"; exit 1; }
done
waitlean; step P2-deps env LEAN_PATH=$LP2 $LEANBIN -DautoImplicit=false $W/Deps.lean
waitlean; step P2-compare-pos env LEAN_PATH=$LP2 $LEANBIN --run /content/work/ours/Compare.lean Work.FourierC588M.UniformFourierChallenge Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
waitlean; step P2-compare-neg env LEAN_PATH=$LP2 $LEANBIN --run /content/work/ours/Compare.lean Work.FourierC588MNeg.UniformFourierChallenge Work.FourierC588M.Main OAI.PowerSaving.transform_mainC588M OAI.PowerSaving.convolution_mainC588M
echo "P2 negative rc=$? (must be nonzero)"
same=0; dif=0; for m in $MODS; do f=$(echo $m | tr . /).olean; if cmp -s $O2/$f $MB/.lake/build/lib/lean/$f; then same=$((same+1)); else dif=$((dif+1)); fi; done
echo "P2 regular files in farm after build: $(find $O2 -type f | wc -l) (expected 140)"
echo "P2 oleans byte-identical to m-build's: $same, differing: $dif"
echo "END $(date -u)"
