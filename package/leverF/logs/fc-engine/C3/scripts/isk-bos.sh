#!/bin/bash
# i-skel: offline re-creation of the bigospeedrun/verifier check (verifier c92b9b6, profile "direct", challenge
# UniformFourier233.lean + Race.bound at the exact bound NUM/10^10) for the instance SUF of my tree.
# fsp-build's fspb-bos-prep.sh + fspb-bos-run.sh merged, with my paths and parameters; nothing else changed:
#   prep: packaging control (Shea's C588 packaged by my script = nmamano's submission, modulo the adapter's comment),
#         stage = proof closure of Work.Fourier<SUF>.Main under Submission/, adapter Submission.lean, VerifierChallenge,
#         VerifierAudit, lakefile / manifest / toolchain from the verifier's profile, Mathlib packages copied;
#   run:  the three phases of verifier.py verify(): lean (comparator, Lean kernel), audit (#print axioms Race.bound),
#         nanoda (comparator with nanoda), each in prep's jail kit.
# Mode "neg": the same stage at bound (NUM+1)/10^10 (one grid step too high), reusing the oleans of stage <X> if it
# exists; only the lean phase runs; expected: the adapter does not compile and the comparator rejects.
# Writes $R/BOS-<X>-OK when every verdict is the expected one.
# usage: isk-bos.sh <X> <SUF> <NUM> [threads, default 3] [pos|neg, default pos] [adapter comment file]
set -uo pipefail
. ${ISK_SCRIPTS:-/content/work/i-skel/scripts}/isk-lib.sh
X=$1; SUF=$2; NUM=$3; THR=${4:-3}; MODE=${5:-pos}; CMT=${6:-}
B=$R/bos-$X; ST=$B/stage; L=$LOGROOT/bos-$X; mkdir -p $L
rm -f $R/BOS-$X-OK
ROOT=Work.Fourier$SUF.Main
if [ $MODE = neg ]; then BOUND=$((NUM+1))/10000000000; else BOUND=$NUM/10000000000; fi
echo "BOS $X: $SUF $ROOT bound $BOUND mode $MODE threads $THR START $(date -u)"
[ "$(git -C $V rev-parse HEAD)" = c92b9b6666ccf5b74c45bd7e664dd3088992f870 ] || fail "verifier clone not at c92b9b6"
echo "verifier clone at c92b9b6; tracked changes: $(git -C $V status --porcelain | wc -l)"
(cd $J && sha256sum -c --quiet HASHES) || fail "jail kit hashes"
rm -rf $B; mkdir -p $ST $B/control/bin $B/scratch
# --- prep -------------------------------------------------------------------------------------------------------------
python3 -I $S/isk_submission.py /content/work/shea/lean/complex588 Work.FourierC588.Main C588 2466337/2000000000 $B/control-c588 || fail "control packaging"
diff -r $B/control-c588/Submission $V/submissions/dft/nmamano-shea-c588/Submission || fail "control packaging differs"
echo "CONTROL PASS: packaging of Shea's C588 identical to nmamano's submission ($(find $B/control-c588/Submission -name '*.lean' | wc -l) files)"
diff <(sed -n '/^import Submission/,$p' $B/control-c588/Submission.lean) <(sed -n '/^import Submission/,$p' $V/submissions/dft/nmamano-shea-c588/Submission.lean) || fail "control adapter body differs"
echo "CONTROL PASS: adapter body identical"
rm -rf $B/control-c588
if [ -n "$CMT" ]; then python3 -I $S/isk_submission.py $T $ROOT $SUF $BOUND $ST $CMT || fail packaging
else python3 -I $S/isk_submission.py $T $ROOT $SUF $BOUND $ST || fail packaging; fi
echo "staged lean files: $(find $ST -name '*.lean' | wc -l); bytes: $(du -sb $ST | cut -f1); largest: $(find $ST -name '*.lean' -printf '%s %p\n' | sort -n | tail -1)"
P=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_file'])")
H=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_sha256'])")
echo "profile challenge: $P  expected sha256 $H"
echo "$H  $V/$P" | sha256sum -c - || fail "profile challenge hash"
python3 -I - "$V/$P" "$ST/VerifierChallenge.lean" "$BOUND" <<'PY'
import sys
from fractions import Fraction
base=open(sys.argv[1],encoding='utf-8').read()
n,d=Fraction(sys.argv[3]).as_integer_ratio()
target=f'''open OAI.PowerSaving Filter Asymptotics

theorem Race.bound : ∃ order solve W, DFTProgram order solve W ∧
    (fun n : ℕ => (W n : ℝ)) =O[atTop]
      (fun n : ℕ => (n : ℝ) * (Real.log (n : ℝ)) ^ (1 - ({n} : ℝ) / {d})) := by
  sorry
'''
open(sys.argv[2],'w',encoding='utf-8').write(base+'\n-- Big O Speedrun exact-bound target; generated from validated integers.\n'+target)
print("challenge bound",n,d)
PY
sha256sum $ST/VerifierChallenge.lean
printf 'import Submission\n#print axioms Race.bound\n' > $ST/VerifierAudit.lean
python3 -I - "$V/profiles/lake-manifest.json" "$ST/lakefile.toml" <<'PY'
import json,sys
result = '''name = "Verifier"
version = "0.1.0"
[leanOptions]
autoImplicit = false
[[require]]
name = "mathlib"
git = "https://github.com/leanprover-community/mathlib4.git"
rev = "d13f23b723b8a846827a245b89c10fc7d3f11612"
'''
m=json.load(open(sys.argv[1]))
mathlib=next(p for p in m['packages'] if p['name']=='mathlib')
result=result.replace('d13f23b723b8a846827a245b89c10fc7d3f11612',mathlib['rev'])
for name in ['Submission','VerifierChallenge','VerifierAudit']:
    result += f'\n[[lean_lib]]\nname = {json.dumps(name)}\n'
open(sys.argv[2],'w').write(result)
PY
echo "leanprover/lean4:v4.34.1" > $ST/lean-toolchain
cp $V/profiles/lake-manifest.json $ST/lake-manifest.json
mkdir -p $ST/.lake/build $ST/.lake/config
cp -a /content/work/shea/lean/complex588/.lake/packages $ST/.lake/packages
if [ $MODE = neg ] && [ -d $R/bos-pos-oleans-$SUF ]; then
  cp -a $R/bos-pos-oleans-$SUF/. $ST/.lake/build/
  find $ST/.lake/build -maxdepth 4 \( -name 'Submission.*' -o -name 'VerifierChallenge.*' -o -name 'VerifierAudit.*' \) -print -delete
  echo "reused the positive stage's Submission/* oleans ($(find $ST/.lake/build -name '*.olean' | wc -l))"
fi
for k in lean nanoda; do
python3 -I - "$B/control/$k.json" "$k" <<'PY'
import json,sys
k=sys.argv[2]
json.dump({'challenge_module':'VerifierChallenge','solution_module':'Submission','theorem_names':['Race.bound'],
 'definition_names':[],'permitted_axioms':['propext','Quot.sound','Classical.choice'],'enable_nanoda':k=='nanoda'},
 open(sys.argv[1],'w'),indent=2,sort_keys=True)
PY
done
echo $THR > $B/control/threads   # lowered before each phase so that my jobs + other heavy Lean stay <= 3
cat > $B/control/bin/lake <<EOS
#!/bin/sh
export LEAN_NUM_THREADS=\$(cat $B/control/threads 2>/dev/null || echo 1)
exec $TC/lake "\$@"
EOS
chmod +x $B/control/bin/lake
chmod 755 $B $B/control $B/control/bin
chown -R 65534:65534 $ST/.lake $B/scratch
du -sh $ST/.lake/packages; echo "PREP DONE $(date -u)"
# --- run --------------------------------------------------------------------------------------------------------------
PHASES="lean audit nanoda"; [ $MODE = neg ] && PHASES=lean
mkdir -p $R/cc/home
ok=1
for ph in $PHASES; do
  wl 14000
  t=$(( MAXHEAVY - $(heavy) )); [ $t -le $THR ] || t=$THR; [ $t -ge 1 ] || t=1; echo $t > $B/control/threads; echo "   bos Lean threads for $ph: $t"
  if [ $ph != audit ]; then find $ST/.lake/build/lib/lean $ST/.lake/build/ir -maxdepth 1 -name 'VerifierChallenge.*' -print -delete 2>/dev/null; fi
  case $ph in
    lean|nanoda) CMD="exec lake env /bin/sh -c 'PATH=$B/control/bin:\$PATH; export PATH; exec comparator $B/control/$ph.json'";;
    audit) CMD="exec lake env lean VerifierAudit.lean";;
  esac
  s=$(date +%s)
  echo "=== $ph START $(date -u) olean files: $(find $ST/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
  echo "cmd: $CMD"
  env -i PATH=$TC:/content/work/cc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=$R/cc/home LANG=C.UTF-8 LEAN_NUM_THREADS=1 JAIL_TMP=$B/scratch \
    /usr/bin/time -v -o $L/$ph.time $J/jail2.sh $ST /bin/sh -c "$CMD" > $L/$ph.log 2>&1
  rc=$?
  echo "=== $ph END rc=$rc wall=$(( $(date +%s)-s ))s $(date -u) maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' $L/$ph.time)"
  tail -4 $L/$ph.log | cut -c1-400
  echo "olean files after: $(find $ST/.lake/build -name '*.olean' 2>/dev/null | wc -l)"
  case $MODE/$ph in
    pos/lean|pos/nanoda) [ $rc = 0 ] && grep -q "Your solution is okay!" $L/$ph.log && echo "VERDICT bos $ph: accepted" || { echo "VERDICT bos $ph: NOT accepted"; ok=0; } ;;
    pos/audit) grep -q "'Race.bound' depends on axioms: \[propext, Classical.choice, Quot.sound\]" $L/$ph.log && echo "VERDICT bos audit: standard axioms" || { echo "VERDICT bos audit: NOT standard"; ok=0; } ;;
    neg/lean) [ $rc != 0 ] && ! grep -q "Your solution is okay!" $L/$ph.log && echo "VERDICT bos negative (bound one grid step too high): rejected (expected)" || { echo "VERDICT bos negative: NOT rejected"; ok=0; } ;;
  esac
done
if [ $MODE = pos ] && [ $ok = 1 ]; then rm -rf $R/bos-pos-oleans-$SUF; mkdir -p $R/bos-pos-oleans-$SUF; cp -a $ST/.lake/build/. $R/bos-pos-oleans-$SUF/; fi
echo "BOS RUN END $(date -u)"
[ $ok = 1 ] && { touch $R/BOS-$X-OK; echo "BOS $X ALL AS EXPECTED"; echo EXIT 0; } || { echo "BOS $X: SOMETHING UNEXPECTED"; echo EXIT 1; exit 1; }
