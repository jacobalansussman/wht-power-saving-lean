#!/bin/bash
# prep (wht45/lever, 19 Oct 2026): finish a full-mode ours-setup.sh run that stopped only at its over-strict sorry check.
# ours-setup.sh grepped ALL 257 overlay files for `sorry` and failed on the 20 hits, which are the two `sorry` lines
# (327, 331) of each UniformFourierChallenge.lean: OpenAI's own statement file (its `sorry` lines are at 327 and 331 too),
# the challenge that Compare.lean checks the solution against.  fsp-build printed the same lines on 14 Oct.
# This script re-checks EVERY result of that run from its logs, with the corrected sorry check (no `sorry` outside the
# challenge files), and only then writes OURS-READY.  No Lean is run.
set -uo pipefail
W=/content/work; L=$W/logs/prep/ours; O=$W/logs/prep/ours.out; T=$W/ours/T/lean/complex588
REF=$W/ref/UniformFourier-fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb.lean
fail() { echo "FINISH FAILED: $*"; exit 1; }
[ -e $W/OURS-READY ] && fail "OURS-READY already exists"
grep -q "mode=full" $O || fail "not a full-mode run"
grep -q "^OURS FAILED: sorry/admit found" $O || fail "the run did not stop at the sorry check"
grep -q "^IN-TREE-OK: 257 overlay files" $O || fail "overlay"
grep -q "^=== build END rc=0" $O || fail "build"
grep -q "^=== nobuild END rc=0" $O && grep -q "All targets up-to-date" $L/nobuild.log || fail "nobuild"
[ "$(grep -c 'depends on axioms: \[propext, Classical.choice, Quot.sound\]' $L/axioms.log)" = 3 ] || fail "final axioms"
[ "$(grep -c 'axioms' $L/axioms-glue.log)" = 16 ] || fail "glue axioms count"
[ -z "$(grep -vE 'axioms: \[(propext|Classical.choice|Quot.sound)(, (propext|Classical.choice|Quot.sound))*\]' $L/axioms-glue.log)" ] || fail "glue axioms"
grep -q "^hostSchedK negative rc=0" $O || fail "hostSchedK negative"
grep -q "^challenge FourierU8X0HQ (1 - 16256629/10^10): byte-identical to my derivation" $O || fail "challenge"
grep -q "^challenge FourierU8X0HQNeg (1 - 16256630/10^10): byte-identical to my derivation" $O || fail "negative challenge"
grep -q "RESULT: PASS" $L/compare-pos.log && grep -q "RESULT: FAIL" $L/compare-neg.log || fail "Compare.lean"
[ ! -s $L/admit-grep.txt ] || fail "admit"
[ "$(grep -n sorry $REF | cut -d: -f1 | tr '\n' ' ')" = "327 331 " ] || fail "reference sorry lines"
grep -vE '/UniformFourierChallenge\.lean:(327|331):  sorry$' $L/sorry-grep.txt > $L/sorry-grep-nonchallenge.txt
[ ! -s $L/sorry-grep-nonchallenge.txt ] || { cat $L/sorry-grep-nonchallenge.txt; fail "sorry outside the challenge files"; }
grep -vE ':[0-9]+:No `sorry`, no `native_decide`\.$' $L/forbidden-grep.txt > $L/forbidden-noncomment.txt
[ ! -s $L/forbidden-noncomment.txt ] || { cat $L/forbidden-noncomment.txt; fail "forbidden construct"; }
[ "$(stat -c %U $T/.lake)" = root ] || fail ".lake owner"
[ -e $W/BASE-READY ] || fail "no BASE-READY"
echo "all checks of the run hold: sorry lines only in the $(cut -d: -f1 $L/sorry-grep.txt | sort -u | wc -l) challenge files (lines 327, 331 = OpenAI's statement), 5 forbidden-grep hits are the comment 'No sorry, no native_decide'"
RES="build OK, --no-build OK, axioms OK, glue axioms printed, hostSchedK negative proved, Compare.lean U8X0HQ positive PASS / negative FAIL, challenges byte-identical to OpenAI's statement with 16256629 / 16256630 (finished by ours-finish.sh)"
{ echo "OURS-READY $(date -u) host=$(hostname) mode=full"; echo "tree /content/work/ours/T (Lake project /content/work/ours/T/lean/complex588)"; echo "$RES"; echo "BASE-READY: $(cat $W/BASE-READY)"; } > $W/OURS-READY
cp $W/ours/pkg/ours-u8x0hq-overlay/meta/OURS-README.txt $W/OURS-README.txt
echo "OURS-READY written"; cat $W/OURS-READY
