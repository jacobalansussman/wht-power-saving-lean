#!/bin/bash
# m-build (13 Oct 2026): stage C588M for an offline re-creation of the bigospeedrun/verifier check, following
# c-base's cb-bos-prep.sh (verifier commit c92b9b6, scripts/verifier.py prepare()) with one change: the Submission
# is not taken from the verifier's repository but packaged from our built tree by mb_submission.py, exactly as
# nmamano packaged Shea's C588 (control below: the same script on Shea's C588 reproduces nmamano's 350 files).
# usage: mb-bos-prep.sh [tree lean/complex588] [bound n/d]
set -eu
T=${1:-/content/work/mb/shea/lean/complex588}
BOUND=${2:-6776303/5000000000}
V=/content/work/bos-verifier
B=/content/work/mb/bos
S=$B/stage
[ "$(git -C $V rev-parse HEAD)" = c92b9b6666ccf5b74c45bd7e664dd3088992f870 ]
rm -rf $B; mkdir -p $S $B/control/bin $B/scratch
# control: our packaging script on Shea's C588 gives nmamano's submission byte for byte (modulo the adapter's comment)
python3 -I /content/work/mb/mb_submission.py /content/work/shea/lean/complex588 Work.FourierC588.Main C588 2466337/2000000000 $B/control-c588
diff -r $B/control-c588/Submission $V/submissions/dft/nmamano-shea-c588/Submission
echo "CONTROL PASS: packaging of Shea's C588 identical to nmamano's submission ($(find $B/control-c588/Submission -name '*.lean' | wc -l) files)"
diff <(sed -n '/^import Submission/,$p' $B/control-c588/Submission.lean) <(sed -n '/^import Submission/,$p' $V/submissions/dft/nmamano-shea-c588/Submission.lean)
echo "CONTROL PASS: adapter body identical"
rm -rf $B/control-c588
python3 -I /content/work/mb/mb_submission.py $T Work.FourierC588M.Main C588M $BOUND $S
echo "staged lean files: $(find $S -name '*.lean' | wc -l); bytes: $(du -sb $S | cut -f1); largest: $(find $S -name '*.lean' -printf '%s %p\n' | sort -n | tail -1)"
P=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_file'])")
H=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_sha256'])")
echo "profile challenge: $P  expected sha256 $H"
echo "$H  $V/$P" | sha256sum -c -
python3 -I - "$V/$P" "$S/VerifierChallenge.lean" "$BOUND" <<'PY'
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
sha256sum $S/VerifierChallenge.lean
printf 'import Submission\n#print axioms Race.bound\n' > $S/VerifierAudit.lean
python3 -I - "$V/profiles/lake-manifest.json" "$S/lakefile.toml" <<'PY'
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
echo "leanprover/lean4:v4.34.1" > $S/lean-toolchain
cp $V/profiles/lake-manifest.json $S/lake-manifest.json
mkdir -p $S/.lake/build $S/.lake/config
cp -a /content/work/shea/lean/complex588/.lake/packages $S/.lake/packages
for k in lean nanoda; do
python3 -I - "$B/control/$k.json" "$k" <<'PY'
import json,sys
k=sys.argv[2]
json.dump({'challenge_module':'VerifierChallenge','solution_module':'Submission','theorem_names':['Race.bound'],
 'definition_names':[],'permitted_axioms':['propext','Quot.sound','Classical.choice'],'enable_nanoda':k=='nanoda'},
 open(sys.argv[1],'w'),indent=2,sort_keys=True)
PY
done
echo 1 > $B/control/threads
cat > $B/control/bin/lake <<EOS
#!/bin/sh
export LEAN_NUM_THREADS=\$(cat $B/control/threads 2>/dev/null || echo 1)
exec /content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lake "\$@"
EOS
chmod +x $B/control/bin/lake
chown -R 65534:65534 $S/.lake $B/scratch
du -sh $S/.lake/packages; echo PREP DONE
