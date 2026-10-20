#!/bin/bash
# c-base: re-create, offline, the core of bigospeedrun/verifier's check of submissions/dft/nmamano-shea-c588
# (verifier commit c92b9b6). Follows scripts/verifier.py prepare(): stage only the Submission .lean files, write the
# verifier's generated challenge (challenge_source: profiles/direct.json "dft" -> challenges/dft/UniformFourier233.lean
# + the Race.bound target), VerifierAudit.lean, its lakefile.toml, lean-toolchain, profiles/lake-manifest.json, and the
# lean/nanoda comparator configs. NOT the verifier's own run: no bwrap/systemd-run/cgroups here; jail.sh instead.
set -eu
V=/content/work/bos-verifier
B=/content/work/bos
S=$B/stage
SUB=submissions/dft/nmamano-shea-c588
[ "$(git -C $V rev-parse HEAD)" = c92b9b6666ccf5b74c45bd7e664dd3088992f870 ]
rm -rf $B; mkdir -p $S $B/control/bin $B/scratch
cat $V/$SUB/manifest.json; echo
# stage: tracked regular .lean files under the root "Submission" (verifier.py prepare(), source_subdir = $SUB)
git -C $V ls-tree -r HEAD $SUB | awk '$1!="100644"{print "NONREGULAR",$0; bad=1} END{exit bad}'
( cd $V/$SUB && git ls-files -z | grep -z '\.lean$' | grep -zE '^Submission(\.lean|/)' | xargs -0 -I{} cp --parents {} $S/ )
echo "staged lean files: $(find $S -name '*.lean' | wc -l)"
# challenge: verifier challenge_source() for {"problem":"dft","bound":"2466337/2000000000"}
P=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_file'])")
H=$(python3 -I -c "import json;print(json.load(open('$V/profiles/direct.json'))['dft']['challenge_sha256'])")
echo "profile challenge: $P  expected sha256 $H"
echo "$H  $V/$P" | sha256sum -c -
python3 -I - "$V/$P" "$S/VerifierChallenge.lean" <<'PY'
import sys
from fractions import Fraction
base=open(sys.argv[1],encoding='utf-8').read()
n,d=Fraction('2466337/2000000000').as_integer_ratio()
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
# packages: same nine revisions as Shea's manifest (checked equal); copy of the cache-fetched packages
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
# lake adapter (stand-in for scripts/lake-shim.c): fixes the thread count, which landrun would otherwise strip
cat > $B/control/bin/lake <<EOS
#!/bin/sh
export LEAN_NUM_THREADS=\$(cat /content/work/bos/control/threads 2>/dev/null || echo 1)
exec /content/work/elan/toolchains/leanprover--lean4---v4.34.1/bin/lake "\$@"
EOS
chmod +x $B/control/bin/lake
cat $S/lakefile.toml; cat $B/control/*.json
chown -R 65534:65534 $S/.lake
du -sh $S/.lake/packages; echo PREP DONE
