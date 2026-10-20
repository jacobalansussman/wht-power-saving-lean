#!/bin/bash
# REGENERATE.sh OUTDIR: make every generated file of the L2979H lever instance (design A: the community's #2979 circuit
# with every head kept) again, from the community's #2979 files (tree/wht45/lever/p-alt/dl/main2979/..., main bb4bccbf),
# the two vendor pins and design A's fill (tree/wht45/allheads/aA-data/in/alloc2979H/FILL.json.gz, 764395fc...), with
# the programs in tree/ (the workspace's own files at the same relative paths; the only changes are the path-only
# patch and one pin in PATCHES.txt), and compare the result byte for byte with ../SHA256SUMS (overlay/...) and the
# controls in ../controls/.  Each generator that writes a figure or a count is also run one unit too high and must
# refuse.  Light Python 3 (standard library) only; a few minutes, under 0.5 GB.
# aA-final (wht45/allheads, 19 Oct 2026): x-final's REGENERATE.sh (lever2979) with steps 3-4 for design A.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); PKG=$(dirname "$HERE")
export GEN_ROOT=$HERE/tree
O=${1:?usage: REGENERATE.sh OUTDIR}; [ -e "$O" ] && { echo "REFUSE: $O exists"; exit 1; }; mkdir -p "$O"; O=$(cd "$O" && pwd)
Q="python3 -I -B"
L=$GEN_ROOT/wht45/lever; X=$GEN_ROOT/wht45/lever2979; BC=$GEN_ROOT/wht43/bound-check
D=$L/p-alt/dl/main2979/research/rank7-joint-endframes; V=$L/g-heads/dl/pr3945/owner/inputs
AH=$GEN_ROOT/wht45/allheads; AT=$AH/aA-data/tools; AL=$AH/aA-data/in/alloc2979H; ME=$AH/aA-final/tools
# aA-data's mk_inst2979H.py reads its data from the folder above its own (out/, rerun/): it runs from a copy of its
# tools/ in $O/aAd (so that nothing is written into tree/), with the regenerated data in $O/aAd/out and $O/aAd/rerun
EXP='{"1": 2445, "2": 706, "3": 232, "4": 108, "5": 18, "8": 8}'
LOG=$O/regenerate.log; : > $LOG
bad=0
run() { echo "+ $*" >> $LOG; "$@" >> $LOG 2>&1 || { echo "FAILED: $*"; bad=1; }; }
refuse() { local want=$1; shift; echo "+ (must exit $want) $*" >> $LOG; "$@" >> $LOG 2>&1; local rc=$?
  if [ $rc = $want ]; then echo "refused as required (exit $rc): ${*: -3}"; else echo "DID NOT REFUSE (exit $rc): $*"; bad=1; fi; }
(cd $HERE/tree && shasum -a 256 -c --quiet ../tree.sha256) || { echo "generators tree differs from tree.sha256"; exit 1; }

echo "1. the two provider words: p-lower's transport and lowering (gauge44 = WORD-base, extra = WORD-extra)"
mkdir -p $O/words
for w in base extra; do
  prov=$([ $w = base ] && echo gauge44 || echo extra)
  run $Q $L/p-lower/tools/transport.py $D/data/WORD-$w.json.gz $O/words/$w.canonical.json.gz $O/words/$w.transport-audit.json
  run $Q $L/p-lower/tools/lower_w.py $O/words/$w.canonical.json.gz $O/words/$w.gcert1w.json.gz $O/words/lower-$w.json \
      --provider $prov --word $D/data/WORD-$w.json.gz --mutants
done
echo "   lowered words: base $(shasum -a 256 $O/words/base.gcert1w.json.gz | cut -c1-16) (x-gate 8f2db1387aefe62b), extra $(shasum -a 256 $O/words/extra.gcert1w.json.gz | cut -c1-16) (x-gate 163466d55a870f90)"
[ "$(shasum -a 256 $O/words/base.gcert1w.json.gz | cut -c1-64)" = 8f2db1387aefe62bb660b664e18a001c996360e205dc2d604a2cd5d307048759 ] && [ "$(shasum -a 256 $O/words/extra.gcert1w.json.gz | cut -c1-64)" = 163466d55a870f90776279b1c9b2540c2d8067139d7b8e326ab4eb2061641d48 ] || { echo "LOWERED WORDS DIFFER"; bad=1; }
echo "2. the word data, glue and price in Lean (p-wgen's generators, instance #2979)"
export KWORD_TOOLS=$BC/k-word/tools
for w in base extra; do N=$([ $w = base ] && echo S2979B || echo S2979E)
  mkdir -p $O/tables-$w $O/price
  run $Q $X/x-words/tools/mk_info.py $O/words/$w.gcert1w.json.gz $O/tables-$w/$w.info.json
  run $Q $X/x-words/tools/pack_ceil2.py $O/tables-$w/$w.info.json $O/tables-$w/$w-mut0-ceil.schedule.json.gz
  run $Q $BC/k-combine/tools/make_sched_x.py $O/words/$w.gcert1w.json.gz $O/tables-$w/$w-mut0-ceil.schedule.json.gz $BC/k-combine/toolkit $O/price/$N
  run $Q $X/x-words/tools/wgen2.py $O/words/$w.gcert1w.json.gz $N $O/w2-$N --expect "$EXP"
  run $Q $X/x-words/tools/gen_global_w.py $O/price/$N.json --toolkit $BC/k-combine/toolkit --out $O/w2-$N --name $N \
      --wgen $O/w2-$N/logs/wgen-$N.json --v2 --glue $X/x-words/tools/glue2.lean.tmpl --expect "$EXP"
  refuse 3 $Q $X/x-words/tools/gen_global_w.py $O/price/$N.json --toolkit $BC/k-combine/toolkit --out $O/never-$N --name $N \
      --wgen $O/w2-$N/logs/wgen-$N.json --v2 --glue $X/x-words/tools/glue2.lean.tmpl --expect "$EXP" --dft 14611616
  run $Q $BC/k-combine/tools/generate_checks_x.py --out $O/w2-$N --original-challenge \
      $BC/k-combine/toolkit/third-party/openai-math/lean/ComparatorChallenges/UniformFourier.lean \
      --comparator $BC/k-combine/toolkit/tools/Compare.lean --name $N
  run $Q $X/x-words/tools/mk_stage.py $N $O/w2-$N $X/x-words/stage/extra none $O/stage
done
run $Q $X/x-words/tools/mk_stage_pair.py $O/stage $X/x-words/stage/pair $O/stage
mkdir -p $O/lean && tar xzf $O/stage/stage-S2979pair.tgz -C $O/lean

echo "3. the host tables (#2979's published fill, for the comparison, and design A's fill), the Lean parts of design A's, the dispatch"
mkdir -p $O/tables
# the residual step reads the lowered words and records the hashes of the lowering reports; those reports name the
# absolute paths of the run that made them, so the step gets the regenerated words (byte-identical, checked) with the
# ORIGINAL reports (tree/.../x-gate/out), after checking that the regenerated reports equal them except for those paths
# and the recorded run time (seconds)
mkdir -p $O/resid-in
for w in base extra; do
  $Q - $O/words/lower-$w.json $X/x-gate/out/lower-$w.json <<'PY2' || { echo "LOWERING REPORT $w DIFFERS"; bad=1; }
import json, sys
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k not in ('canonical', 'word', 'out', 'native', 'seconds')}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
a, b = (json.load(open(f)) for f in sys.argv[1:3])
assert strip(a) == strip(b), 'reports differ beyond the recorded paths'
print('lowering report equal to the original except for the recorded paths and run time:', sys.argv[2].rsplit('/', 1)[1])
PY2
  cp $O/words/$w.gcert1w.json.gz $O/resid-in/; cp $X/x-gate/out/lower-$w.json $O/resid-in/
done
run $Q $X/x-tables/tools/mk_resid2979.py $O/resid-in $O/tables/residuals-2979.json
# Published copy (sixth revision of the repository, 20 Oct 2026). The two lowering reports in tree/.../x-gate/out name
# the machine paths of the run that made them; those paths were shortened for publication, which changed the two
# files' sha256 (tree.sha256 gives the published bytes). The residual file records the sha256 of each report, so the
# two recorded hashes are set back to those of the reports as the run wrote them (f219f479..., 5e5d7945...); the pins
# after this step then check everything exactly as before. No generator is changed.
$Q - $O/tables/residuals-2979.json <<'PY3' || { echo "RESIDUAL REPORT HASHES NOT RESTORED"; bad=1; }
import json, sys
p = sys.argv[1]
d = json.load(open(p))
assert set(d['source']) == {'base.gcert1w.json.gz', 'lower-base.json', 'extra.gcert1w.json.gz', 'lower-extra.json'}
d['source']['lower-base.json'] = 'f219f479f8b72d772648715709686c171744580bb35bc23757bc99b6cac10a2d'
d['source']['lower-extra.json'] = '5e5d794595c9d370c9212900bdb915d002b88833de9ea41d0f29982010924012'
json.dump(d, open(p, 'w'))
print('   residuals: the two lowering reports recorded with their sha256 as the run wrote them')
PY3
run $Q $X/x-tables/tools/own_to_host2979.py $D/allocation $D/data $V $O/tables/residuals-2979.json $O/tables/H2979
[ "$(shasum -a 256 $O/tables/residuals-2979.json | cut -c1-64)" = 05503c4534460602ed4181f555924fc7529c0bef7e7fb98a05ba9dc99c1275d5 ] \
  && echo "   residuals = x-tables' (05503c45...)" || { echo "RESIDUALS DIFFER FROM X-TABLES'"; bad=1; }
# design A: the same converter with only the FILL pin changed (aA-data's own_to_host2979H.py), on design A's allocation
# folder (its FILL; PROFILE, ORDER, DISPATCH, PROVIDERS are the community's #2979 files, byte for byte)
for f in PROFILE.json ORDER.json DISPATCH.json PROVIDERS.json; do cmp -s $AL/$f $D/allocation/$f || { echo "design A's $f is not #2979's"; bad=1; }; done
[ "$(shasum -a 256 $AL/FILL.json.gz | cut -c1-64)" = 764395fcc2d883e9dc1ecfedce2131c745dc7ad1b9c41d7b52ef54da16d7f71c ] && echo "   design A's FILL 764395fc..." || { echo "DESIGN A FILL DIFFERS"; bad=1; }
mkdir -p $O/aAd/out $O/aAd/rerun/disp; cp -R $AT $O/aAd/tools
run $Q $O/aAd/tools/own_to_host2979H.py $AL $D/data $V $O/tables/residuals-2979.json $O/aAd/out/H2979H
run $Q $ME/fill_indep2979H.py $O/tables/H2979-hosttab.json.gz $O/aAd/out/H2979H-hosttab.json.gz $O/aAd/out/fill-independence.json
# the head count (16,576 = 8,288 + 8,288 literal heads; 11,614 other rank-48 blocks = 6,342 + 5,272): heads2979.py
# prints ALL OK or SOME CHECK FAILED (exit 0 either way), so the verdict line is checked; one head more must fail
$Q $O/aAd/tools/heads2979.py $O/aAd/out/H2979H-hosttab.json.gz $AL/FILL.json.gz $V/TEMPLATE.json.gz - \
    $O/aAd/out/heads2979H.json.gz $O/aAd/out/heads2979H-summary.json 8288,8288,6342,5272 > $O/tables/heads.out 2>&1
tail -1 $O/tables/heads.out | grep -qx 'ALL OK' && echo "head count 16,576 (8,288 + 8,288): ALL OK" || { echo "HEAD COUNT CHECK FAILED"; bad=1; }
$Q $O/aAd/tools/heads2979.py $O/aAd/out/H2979H-hosttab.json.gz $AL/FILL.json.gz $V/TEMPLATE.json.gz - \
    $O/tables/never-heads.json.gz $O/tables/never-heads-summary.json 8289,8288,6342,5272 > $O/tables/heads-neg.out 2>&1
tail -1 $O/tables/heads-neg.out | grep -qx 'SOME CHECK FAILED' && echo "refused as required: one head more (8,289)" || { echo "ONE HEAD MORE NOT REFUSED"; bad=1; }
run $Q $O/aAd/tools/mk_parts.py $O/aAd/out/H2979H-hosttab.json.gz H2979H 25 $O/lean/Work/XHostT/Data/H2979H
mkdir -p $O/controls-made; mv $O/lean/Work/XHostT/Data/H2979H/NegGlob.lean $O/controls-made/NegGlob.lean   # a negative control: compared with ../controls/
refuse 3 $Q $O/aAd/tools/mk_parts.py $O/aAd/out/H2979H-hosttab.json.gz H2979Hnload 25 $O/never-load --perturb load
refuse 3 $Q $O/aAd/tools/mk_parts.py $O/aAd/out/H2979H-hosttab.json.gz H2979Hnwin 25 $O/never-win --perturb win
mkdir -p $O/tables/disp $O/lean/Work/XHostT/Data/D2979
run $Q $X/x-tables/tools/mk_disp2979.py $D/allocation/DISPATCH.json $D/allocation/ORDER.json 74 416 250 $O/tables/disp/Disp.lean
cp $O/tables/disp/Disp.lean $O/lean/Work/XHostT/Data/D2979/Disp.lean; cp $O/tables/disp/Disp.lean $O/aAd/rerun/disp/Disp.lean
mkdir -p $O/never-disp; refuse 3 $Q $X/x-tables/tools/mk_disp2979.py $D/allocation/DISPATCH.json $D/allocation/ORDER.json 74 417 249 $O/never-disp/Disp.lean

echo "4. the rate fact at the figure (e-r45m1's generator through aA-data's wrapper), the instance modules (aA-data's checked"
echo "   renaming), the statement chain, the warp heads-first files, the engine theorem and the controls"
run $Q $O/aAd/tools/mk_profile2979H.py $O/aAd/out/H2979H-audit.json $O/aAd/out/heads2979H-summary.json $O/aAd/profile2979H.json
cmp -s $O/aAd/profile2979H.json $AT/profile2979H.json && echo "   profile2979H.json regenerated byte-identical" || { echo "PROFILE2979H DIFFERS"; bad=1; }
run $Q $L/e-r45m1/tools/mk_rate_mixed.py scan 2979H=$AT/profile2979H.json
run $Q $O/aAd/tools/mk_rate2979H.py $AT/profile2979H.json 17564707 $O/lean $O/aAd/out/rate-2979H.json
mkdir -p $O/never-rate; refuse 2 $Q $O/aAd/tools/mk_rate2979H.py $AT/profile2979H.json 17564708 $O/never-rate $O/never-rate/never.json
mv $O/lean/Work/Mixed/Rate2979HSharp.lean $O/controls-made/Rate2979HSharp.lean    # a side check (../controls/side/)
run $Q $O/aAd/tools/mk_inst2979H.py $O/lean
mv $O/lean/Work/XHostT/Data/H2979H/NegHeads.lean $O/controls-made/NegHeads.lean   # a negative control: compared with ../controls/
for f in Work/XHostT/Inst/S2979HHB.lean Work/Mixed/R5Joint2979H.lean; do mkdir -p $O/controls-made/side/$(dirname $f); mv $O/lean/$f $O/controls-made/side/$f; done
(cd $AH/aA-data && run $Q tools/mk_l2979H.py final $O/final-L2979H)
(cd $AH/aA-data && refuse 2 $Q tools/mk_l2979H.py final $O/never-L2979H --num 17564708)
cp -R $O/final-L2979H/stage/. $O/lean/
run $Q $X/x-inst/tools/mk_hfw2979.py $L/p-bridge/lean2 $O/lean
run $Q $ME/mk_prov2979H.py $O/lean
run $Q $ME/mk_controls2979H.py $O/controls-made/top $O/controls-made/Rate2979HSharp.lean $O/controls-made/NegHeads.lean $O/controls-made/NegGlob.lean
for f in Heads2979HNeg.lean Rate2979HSharpCtl.lean NegHeads.lean NegGlob.lean; do
  cmp -s $O/controls-made/top/$f $PKG/controls/$f && echo "   control $f regenerated byte-identical" || { echo "CONTROL $f DIFFERS"; bad=1; }; done
for f in Work/Mixed/Rate2979HSharp.lean Work/XHostT/Inst/S2979HHB.lean Work/Mixed/R5Joint2979H.lean; do
  g=$O/controls-made/side/$f; [ -e $g ] || g=$O/controls-made/$(basename $f)
  cmp -s $g $PKG/controls/side/$f && echo "   side check $f regenerated byte-identical" || { echo "SIDE CHECK $f DIFFERS"; bad=1; }; done

echo "5. comparison with ../SHA256SUMS"
$Q - "$PKG/SHA256SUMS" "$O/lean" "$O/compare.tsv" <<'PY' || bad=1
import hashlib, os, sys
sums, root, out = sys.argv[1:4]
want = {}
for line in open(sums):
    h, p = line.rstrip('\n').split('  ', 1)
    if p.startswith('overlay/'):
        want[p[len('overlay/'):]] = h
made = {}
for dp, dn, fn in os.walk(root):
    for f in fn:
        q = os.path.join(dp, f)
        made[os.path.relpath(q, root)] = hashlib.sha256(open(q, 'rb').read()).hexdigest()
same = sorted(p for p in made if p in want and made[p] == want[p])
diff = sorted(p for p in made if p in want and made[p] != want[p])
extra = sorted(p for p in made if p not in want)
notmade = sorted(p for p in want if p not in made)
with open(out, 'w') as f:
    for p in same: f.write('same\t' + p + '\n')
    for p in diff: f.write('DIFFERS\t' + p + '\n')
    for p in extra: f.write('made, not in the overlay (outside the closure)\t' + p + '\n')
    for p in notmade: f.write('in the overlay, not generated here (hand-written, or generic: see README.txt)\t' + p + '\n')
print(f'regenerated overlay files byte-identical: {len(same)}; DIFFERENT: {len(diff)}; made but outside the closure: {len(extra)}; '
      f'overlay files not generated here: {len(notmade)} (list: {out})')
for p in diff: print('DIFFERS', p)
sys.exit(1 if diff or not same else 0)
PY
if [ $bad = 0 ]; then echo "REGENERATE: ALL AS EXPECTED (log $LOG)"; exit 0; else echo "REGENERATE: SOMETHING FAILED (log $LOG)"; exit 1; fi
