#!/bin/bash
# i-skel: the checks of a built Fourier instance SUF at DFT figure 1 - NUM/10^10 (fsp-build fspb-u8x0hq.sh and
# e-r45m1 mxt-build.sh, with parameters): axioms of the chain's Axioms.lean (+ an optional extra axioms file), sorry /
# forbidden-construct greps over the installed files, both challenge files re-derived from OpenAI's statement (NUM and
# NUM+1) and compared byte for byte, Compare.lean positive (must PASS) and negative (must FAIL).
# Writes $R/CHECKS-<SUF>-OK when every verdict is the expected one.
# usage: isk-checks.sh <SUF> <NUM> [<TAG of the install whose files are grepped>] [<extra axioms .lean file>]
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
SUF=$1; NUM=$2; TAG=${3:-}; EXTRA=${4:-}
L=$LOGROOT/$(lower $SUF); mkdir -p $L
rm -f $R/CHECKS-$SUF-OK
cd $T || fail "no $T"
echo "CHECKS $SUF NUM=$NUM START $(date -u)"
ok=1
# 1. axioms
wl; step axioms lake env lean Work/Fourier$SUF/Axioms.lean; cat $L/axioms.log
axioms_ok $L/axioms.log && echo "axioms $SUF: standard only" || { echo "AXIOMS NOT STANDARD"; ok=0; }
if [ -n "$EXTRA" ]; then wl; step axioms-extra lake env lean $EXTRA; cat $L/axioms-extra.log
  axioms_ok $L/axioms-extra.log && ! grep -q "error" $L/axioms-extra.log && echo "extra axioms: standard only" || { echo "EXTRA AXIOMS NOT STANDARD"; ok=0; }; fi
# 2. greps over the installed files (the challenge files carry OpenAI's two `sorry` lines by design: checked by 3.)
if [ -n "$TAG" ] && [ -e $R/in-$TAG.files ]; then FL=$R/in-$TAG.files; else ls Work/Fourier$SUF/*.lean Work/Fourier${SUF}Neg/*.lean > $L/files.txt; FL=$L/files.txt; fi
grep -v 'UniformFourierChallenge.lean$' $FL > $L/checked-lean-files.txt
echo "files grepped: $(wc -l < $L/checked-lean-files.txt) (from $FL)"
(xargs grep -Hnw -E 'sorry|admit' < $L/checked-lean-files.txt | grep -v '`sorry`') > $L/sorry-grep.txt; echo "sorry/admit lines: $(wc -l < $L/sorry-grep.txt)"; cat $L/sorry-grep.txt
# (a construct named inside backquotes in a comment, e.g. "No `sorry`, no `native_decide`.", is not a use: such spans
#  are removed before matching, as the sorry grep drops lines with `sorry`)
FPAT='^\s*(private |protected |noncomputable )?axiom |native_decide|implemented_by|@\[extern|unsafe |ofReduceBool|ofReduceNat|skipKernelTC|debug\.skipKernelTC|Lean\.Elab|elab_rules|macro_rules'
(xargs grep -HnE "$FPAT" < $L/checked-lean-files.txt) > $L/forbidden-grep-raw.txt
FPAT="$FPAT" perl -ne '($a=$_) =~ s/`[^`]*`//g; $a =~ s/^[^:]*:[0-9]+://; print if $a =~ /$ENV{FPAT}/' $L/forbidden-grep-raw.txt > $L/forbidden-grep.txt
echo "forbidden-construct mentions inside backquotes (not uses): $(( $(wc -l < $L/forbidden-grep-raw.txt) - $(wc -l < $L/forbidden-grep.txt) ))"
echo "forbidden-construct lines: $(wc -l < $L/forbidden-grep.txt)"; head $L/forbidden-grep.txt
[ -s $L/sorry-grep.txt ] && ok=0; [ -s $L/forbidden-grep.txt ] && ok=0
# 3. challenge files from OpenAI's statement (k-build / lg-a1 derive_challenges.sh recipe)
echo "OpenAI ref: $REF_SHA $REF"
NEG=$((NUM+1))
for spec in "$NUM Fourier$SUF" "$NEG Fourier${SUF}Neg"; do set -- $spec
  perl -Mutf8 -CSD -pe "s{1 - 1/\(10\^\(13:ℕ\)\)}{1 - $1/(10^(10:ℕ))}; s/\b(decimalExponent|decimalTime|TimeBounds|DFTGoal|ConvGoal|transform_main|convolution_main)\b/\$1$SUF/g" $REF > $L/derived-$2.lean
  if cmp -s $L/derived-$2.lean Work/$2/UniformFourierChallenge.lean; then echo "challenge $2 (1 - $1/10^10): byte-identical to my derivation from OpenAI's statement"
  else echo "CHALLENGE $2 DIFFERS"; ok=0; fi; done
diff <(sed 's/[0-9]\{6,\}/N/g' Work/Fourier$SUF/UniformFourierChallenge.lean) <(sed 's/[0-9]\{6,\}/N/g' Work/Fourier${SUF}Neg/UniformFourierChallenge.lean) > /dev/null \
  && echo "positive and negative challenge differ only in the numerator" || { echo "CHALLENGES DIFFER BEYOND THE NUMERATOR"; ok=0; }
# 4. Compare.lean positive / negative (skipped with ISK_NOCOMPARE=1: the hypothesis-mode chain, whose final theorems take
#    the engine fact as an argument, cannot match the challenge by design)
if [ "${ISK_NOCOMPARE:-0}" = 1 ]; then
  echo "Compare.lean skipped (ISK_NOCOMPARE=1)"
  if [ $ok = 1 ]; then touch $R/CHECKS-$SUF-NOCMP-OK; echo "CHECKS $SUF (no Compare) ALL AS EXPECTED $(date -u)"; echo EXIT 0; exit 0; else echo "CHECKS $SUF: SOMETHING UNEXPECTED"; echo EXIT 1; exit 1; fi
fi
[ "$(sha256sum tools/Compare.lean | cut -c1-64)" = $COMPARE_SHA ] && echo "Compare.lean e48d6d3b (= /content/work/ours/Compare.lean)" || { echo "COMPARE.LEAN CHANGED"; ok=0; }
wl 20000; step compare-pos lake env lean --run tools/Compare.lean Work.Fourier$SUF.UniformFourierChallenge Work.Fourier$SUF.Main OAI.PowerSaving.transform_main$SUF OAI.PowerSaving.convolution_main$SUF; echo "compare-pos rc=$?"; grep -E "RESULT" $L/compare-pos.log
wl 20000; step compare-neg lake env lean --run tools/Compare.lean Work.Fourier${SUF}Neg.UniformFourierChallenge Work.Fourier$SUF.Main OAI.PowerSaving.transform_main$SUF OAI.PowerSaving.convolution_main$SUF; echo "compare-neg rc=$? (must be nonzero)"; grep -E "RESULT|differs" $L/compare-neg.log | head -5
grep -q "RESULT: PASS" $L/compare-pos.log && echo "VERDICT Compare positive: PASS" || { echo "VERDICT Compare positive: NOT PASS"; ok=0; }
grep -q "RESULT: FAIL" $L/compare-neg.log && ! grep -q "RESULT: PASS" $L/compare-neg.log && echo "VERDICT Compare negative: FAIL (expected)" || { echo "VERDICT Compare negative: NOT FAIL"; ok=0; }
if [ $ok = 1 ]; then touch $R/CHECKS-$SUF-OK; echo "CHECKS $SUF ALL AS EXPECTED $(date -u)"; echo EXIT 0; else echo "CHECKS $SUF: SOMETHING UNEXPECTED $(date -u)"; echo EXIT 1; exit 1; fi
