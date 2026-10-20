#!/bin/bash
# prep (wht45/lever, 19 Oct 2026): OUR best Lean stack GXD.U8X0HQ (DFT of every length, power saving 1 - 16256629/10^10,
# i.e. 1.6256629e-3) rebuilt on top of the prepared base, as fsp-build did on 14 Oct (REPRODUCE.sh, fspb-control.sh,
# fspb-u8x0hq.sh), in ONE step instead of the old chain of copies:
#   /content/work/ours/T = cp -a /content/work/shea (Ryan Shea's tree at 5239484, built by lean-base-setup.sh)
#                          + the overlay ours-u8x0hq-overlay.tgz: 257 source files (fsp-build's 234 base sources =
#                          k-combine U8XM overlay + fs-build u8xmhh-src + xhost-src + fsp-cc stage-u8x0 + stage-u8x0hh,
#                          then fsp-build's 23 stage-u8x0hq files) and lakefile.lean + `lean_lib XExit` + `lean_lib CSX`.
# The copy's .lake stays root-owned; the base tree is only read.  Lake as root in MY tree only.
# usage: ours-setup.sh full|quick   (full adds: Compare.lean positive PASS / negative FAIL, challenge derivation,
#                                    hostSchedK negative, glue axioms, source greps; quick: build + --no-build + axioms)
# log: /content/work/logs/prep/ours.out (steps in /content/work/logs/prep/ours/)
set -uo pipefail
. /content/work/env.sh
W=/content/work; O=$W/ours; L=$W/logs/prep/ours; mkdir -p $L
MODE=${1:-quick}
PKG=$O/pkg; TT=$O/T; T=$TT/lean/complex588
TGZ=$O/in/ours-u8x0hq-overlay.tgz
REF=$W/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean
heavy() { ps -C lean -o rss= | awk '$1>2000000' | wc -l; }
avail() { free -m | awk '/Mem:/ {print $7}'; }
wl() { local need=${1:-14000} w=0; while [ "$(heavy)" -ge 3 ] || [ "$(avail)" -lt $need ]; do sleep 10; w=$((w+10)); done; echo "   (slot: waited ${w}s; heavy lean now $(heavy), avail $(avail) MB)"; }
step() { local lab=$1; shift; local s=$(date +%s); echo "=== $lab START $(date -u +%T)"
  /usr/bin/time -v -o $L/$lab.time "$@" > $L/$lab.log 2>&1; local rc=$?
  echo "=== $lab END rc=$rc wall=$(( $(date +%s)-s ))s maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$lab.time) $(date -u +%T)"; return $rc; }
fail() { echo "OURS FAILED: $*"; echo "OURS END $(date -u)"; echo "EXIT 1"; exit 1; }
echo "OURS START $(date -u) mode=$MODE host=$(hostname)"
# 0. the base: wait until Shea's tree is BUILT (cb-build.sh printed "BUILD rc=0"; afterwards cb-build.sh only reads the
#    tree: axioms, Compare runs, --no-build).  The copy starts then; OURS-READY is written only after the base controls
#    passed too (BASE-READY's three conditions) and BASE-READY exists (comparator tools).
BL=$W/logs/base-build.log
built() { [ -e $W/BASE-READY ] || grep -q "^BUILD rc=0" $BL 2>/dev/null; }
buildfail() { grep -E "^BUILD rc=" $BL 2>/dev/null | grep -qv "rc=0"; }
baseok() { [ -e $W/BASE-READY ] && return 0
  grep -q "^BUILD END" $BL 2>/dev/null && grep -q "BUILD rc=0" $BL && grep -q "STEP ours-compare-fourier END rc=0.*RESULT: PASS" $BL && grep -q "STEP ours-compare-negative END rc=1.*RESULT: FAIL" $BL; }
basefail() { grep -q "^== END" $W/logs/base-all.log 2>/dev/null && [ ! -e $W/BASE-READY ]; }
w=0; until built; do { buildfail || basefail; } && fail "base build failed"; sleep 15; w=$((w+15)); done
echo "base tree built (waited ${w}s; controls $(baseok && echo passed || echo 'still running'))"
[ -e $TT ] && fail "$TT exists"
# 1. overlay package: tarball hash, per-file hashes, no extra files
echo "fc787b0228792ad3f060056956a2faf34e64e863b181f0fe147575bb15bca8da  $TGZ" | sha256sum -c - || fail "tarball hash"
rm -rf $PKG; mkdir -p $PKG && tar xzf $TGZ -C $PKG --no-same-owner || fail "untar"
P=$PKG/ours-u8x0hq-overlay
(cd $P/tree && sha256sum -c --quiet ../SHA256SUMS.tree) || fail "overlay file hashes"
n=$(cd $P/tree && find . -type f | wc -l); [ "$n" = "$(wc -l < $P/SHA256SUMS.tree)" ] || fail "tree has $n files, list has $(wc -l < $P/SHA256SUMS.tree)"
echo "OVERLAY-OK ($n files)"
# 2. my copy of the base tree
s=$(date +%s); cp -a $W/shea $TT || fail "copy"; echo "copy done in $(( $(date +%s)-s ))s; $(du -sh $TT | cut -f1)"
cd $T || fail "no $T"
[ "$(stat -c %U .lake)" = root ] && [ "$(stat -c %U .lake/packages/mathlib)" = root ] || fail ".lake not root-owned"
[ "$(git -C $TT rev-parse HEAD)" = 52394849c7ecf3412d39b925d75fc7fdb98994d8 ] || fail "base HEAD $(git -C $TT rev-parse HEAD)"
[ -z "$(git -C $TT status --porcelain --untracked-files=no)" ] || fail "base has tracked changes"
echo "8927b3ee27366af0a026fbaefbc0c5d87fcbfc9138a68ab0938bf71e6cc9f8b5  lakefile.lean" | sha256sum -c --quiet - || fail "base lakefile differs"
echo "mathlib: $(git -C .lake/packages/mathlib rev-parse HEAD); packages: $(ls .lake/packages | tr '\n' ' ')"
(cd $P/tree && find . -type f) | while read f; do if [ -e "$T/$f" ]; then echo "REFUSE existing $f"; exit 1; fi; done || fail "overlay would overwrite a base file"
cp -a $P/tree/. $T/ && cp $P/meta/lakefile.lean.ours $T/lakefile.lean || fail "overlay copy"
chown -R root:root $T/Work $T/CSX $T/XExit $T/comparator $T/lakefile.lean
echo "b8d12a30852ec8720ce77f37295f531e92a68062f56e08ee0eacbb329baab866  lakefile.lean" | sha256sum -c --quiet - || fail "lakefile"
sha256sum -c --quiet $P/SHA256SUMS.tree || fail "in-tree hashes"
echo "IN-TREE-OK: $(wc -l < $P/SHA256SUMS.tree) overlay files + lakefile ($(tail -2 lakefile.lean | tr '\n' ' '))"
echo "tracked changes: $(git -C $TT status --porcelain --untracked-files=no | tr '\n' ' ')"
# 3. Lake build of the U8X0HQ targets (and the hosted network modules HBank / HEnd, which other tracks import)
TG="Work.BlockAccounting.Moment4 Work.XHost.U8X0HQ Work.FourierU8X0HQ.Main Work.FourierU8X0HQ.UniformFourierChallenge Work.FourierU8X0HQNeg.UniformFourierChallenge Work.FourierU8X0HQ.Axioms Work.XHost.HBank Work.XHost.HEnd"
wl; step build env LEAN_NUM_THREADS=2 lake build $TG; rc=$?
grep -E "✖|error:|Build completed" $L/build.log | head -20
grep -E "Built .* \([0-9.]+(s|m)\)" $L/build.log | sed -E 's/.*Built ([^ ]+) \(([0-9.]+)(s|m)\).*/   \2\3 \1/' | sort -g -r | head -10
[ $rc -eq 0 ] || fail "lake build rc=$rc"
wl 8000; step nobuild lake build --no-build $TG; rc=$?; tail -1 $L/nobuild.log; [ $rc -eq 0 ] || fail "not up to date"
wl; step axioms lake env lean Work/FourierU8X0HQ/Axioms.lean; cat $L/axioms.log
[ "$(grep -c "depends on axioms: \[propext, Classical.choice, Quot.sound\]" $L/axioms.log)" = 3 ] || fail "axioms"
RES="build OK, --no-build OK, axioms OK"
if [ "$MODE" = full ]; then
  wl; step axioms-glue lake env lean $P/meta/AxU8X0HQ.lean; cat $L/axioms-glue.log
  wl; step hsneg lake env lean $P/meta/HsU8X0HQNeg.lean; echo "hostSchedK negative rc=$? (0 = the kernel proved '= false')"; tail -3 $L/hsneg.log
  (cd $P/tree && find Work CSX XExit -name '*.lean' | sort -u) > $L/checked-lean-files.txt
  (xargs grep -nw sorry < $L/checked-lean-files.txt | grep -v '`sorry`') > $L/sorry-grep.txt; echo "sorry lines: $(wc -l < $L/sorry-grep.txt)"
  (xargs grep -nwE 'admit' < $L/checked-lean-files.txt) > $L/admit-grep.txt; echo "admit lines: $(wc -l < $L/admit-grep.txt)"
  (xargs grep -nE '^\s*(private |protected |noncomputable )?axiom |native_decide|implemented_by|@\[extern|unsafe |ofReduceBool|ofReduceNat|skipKernelTC|debug\.skipKernelTC|Lean\.Elab|elab_rules|macro_rules' < $L/checked-lean-files.txt) > $L/forbidden-grep.txt; echo "forbidden-construct lines: $(wc -l < $L/forbidden-grep.txt)"; head $L/forbidden-grep.txt
  echo "OpenAI ref: $(sha256sum $REF)"
  for spec in "16256629 FourierU8X0HQ" "16256630 FourierU8X0HQNeg"; do set -- $spec
    bash $P/meta/derive_challenges.sh $REF $1 U8X0HQ > $L/derived-$2.lean
    cmp -s $L/derived-$2.lean Work/$2/UniformFourierChallenge.lean && echo "challenge $2 (1 - $1/10^10): byte-identical to my derivation" || fail "CHALLENGE $2 DIFFERS"; done
  echo "Compare.lean: $(sha256sum tools/Compare.lean | cut -c1-64) (ours: $(sha256sum $W/ours/Compare.lean | cut -c1-64))"
  wl 20000; step compare-pos lake env lean --run tools/Compare.lean Work.FourierU8X0HQ.UniformFourierChallenge Work.FourierU8X0HQ.Main OAI.PowerSaving.transform_mainU8X0HQ OAI.PowerSaving.convolution_mainU8X0HQ; echo "compare-pos rc=$?"; grep -E "RESULT" $L/compare-pos.log
  wl 20000; step compare-neg lake env lean --run tools/Compare.lean Work.FourierU8X0HQNeg.UniformFourierChallenge Work.FourierU8X0HQ.Main OAI.PowerSaving.transform_mainU8X0HQ OAI.PowerSaving.convolution_mainU8X0HQ; echo "compare-neg rc=$? (must be nonzero)"; grep -E "RESULT|differs" $L/compare-neg.log | head -5
  grep -q "RESULT: PASS" $L/compare-pos.log && grep -q "RESULT: FAIL" $L/compare-neg.log || fail "Compare.lean controls"
  # the only allowed `sorry` lines: 327 and 331 of each UniformFourierChallenge.lean (OpenAI's statement file, the
  # challenge Compare.lean checks against; its own sorry lines are 327 and 331)
  grep -vE '/UniformFourierChallenge\.lean:(327|331):  sorry$' $L/sorry-grep.txt > $L/sorry-grep-nonchallenge.txt
  { [ -s $L/sorry-grep-nonchallenge.txt ] || [ -s $L/admit-grep.txt ]; } && fail "sorry outside the challenge files, or admit"
  RES="$RES, glue axioms printed, Compare.lean U8X0HQ positive PASS / negative FAIL, challenges byte-identical to OpenAI's statement with 16256629 / 16256630"
fi
[ "$(stat -c %U .lake)" = root ] || fail ".lake owner changed"
# 4. OURS-READY only after the base controls passed and BASE-READY (comparator tools) exists
w=0; until baseok; do basefail && fail "base controls did not pass"; sleep 20; w=$((w+20)); done; echo "base controls passed (waited ${w}s)"
w=0; until [ -e $W/BASE-READY ]; do basefail && { echo "WARNING: comparator tools step failed; OURS-READY written without BASE-READY"; break; }; sleep 20; w=$((w+20)); done
echo "BASE-READY waited ${w}s"
{ echo "OURS-READY $(date -u) host=$(hostname) mode=$MODE"; echo "tree /content/work/ours/T (Lake project /content/work/ours/T/lean/complex588)"; echo "$RES"; echo "BASE-READY: $(cat $W/BASE-READY 2>/dev/null || echo MISSING)"; } > $W/OURS-READY
cp $P/meta/OURS-README.txt $W/OURS-README.txt
echo "OURS-READY written"; cat $W/OURS-READY
echo "OURS END $(date -u)"
echo "EXIT 0"
