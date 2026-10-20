#!/bin/bash
# REGENERATE.sh OUTDIR: make every generated file of the L2979Fc lever instance (A-fc: the community's #2979 words with
# the re-solved fill u1-2979sep, its 16,576 literal heads cut by all 48 axes, nothing split) again, from the community's
# #2979 files (tree/wht45/lever/p-alt/dl/main2979/..., main bb4bccbf), the two vendor pins and A-fc's fill
# (tree/wht45/afc/afc-data/in/alloc2979Fc/FILL.json.gz, 7fc766c6..., the same file as ../fill/FILL.json.gz), with the
# programs in tree/ (the workspace's own files at the same relative paths; the only changes are the path-only patch and
# the pins in PATCHES.txt), and compare the result byte for byte with ../SHA256SUMS (overlay/...) and the controls in
# ../controls/.  Each generator that writes a figure, a count or a table is also run one step too far and must refuse.
# Light Python 3 (standard library) only; a few minutes, under 0.5 GB.
# afc-final (wht45/afc, 20 Oct 2026): design A's REGENERATE.sh (aA-final, package/leverA) with steps 3-5 for A-fc:
# the words, the word stage, the residuals and the dispatch are design A's and must come out the same; the host table,
# the instance modules, the rate fact, the statement chain, the providers and the controls are A-fc's (afc-data's,
# afc-cut's and afc-final's generators).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); PKG=$(dirname "$HERE")
export GEN_ROOT=$HERE/tree
O=${1:?usage: REGENERATE.sh OUTDIR}; [ -e "$O" ] && { echo "REFUSE: $O exists"; exit 1; }; mkdir -p "$O"; O=$(cd "$O" && pwd)
Q="python3 -I -B"
L=$GEN_ROOT/wht45/lever; X=$GEN_ROOT/wht45/lever2979; BC=$GEN_ROOT/wht43/bound-check
D=$L/p-alt/dl/main2979/research/rank7-joint-endframes; V=$L/g-heads/dl/pr3945/owner/inputs
AH=$GEN_ROOT/wht45/allheads; ME=$AH/aA-final/tools
AF=$GEN_ROOT/wht45/afc; AT=$AF/afc-data/tools; AL=$AF/afc-data/in/alloc2979Fc; CT=$AF/afc-cut/tools; FT=$AF/afc-final/tools
# afc-data's mk_inst2979Fc.py reads its data from the folder above its own (in/, out/, rerun/): it runs from a copy of
# its tools/ in $O/afd (so that nothing is written into tree/), with the regenerated data in $O/afd/out and $O/afd/rerun
EXP='{"1": 2445, "2": 706, "3": 232, "4": 108, "5": 18, "8": 8}'
LOG=$O/regenerate.log; : > $LOG
bad=0
run() { echo "+ $*" >> $LOG; "$@" >> $LOG 2>&1 || { echo "FAILED: $*"; bad=1; }; }
refuse() { local want=$1; shift; echo "+ (must exit $want) $*" >> $LOG; "$@" >> $LOG 2>&1; local rc=$?
  if [ $rc = $want ]; then echo "refused as required (exit $rc): ${*: -3}"; else echo "DID NOT REFUSE (exit $rc): $*"; bad=1; fi; }
shaf() { shasum -a 256 "$1" | cut -c1-64; }
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

echo "3. the host tables (#2979's published fill, for the comparison, and A-fc's fill), the Lean parts of A-fc's, the dispatch"
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
# A-fc: the same converter with only the FILL and PROFILE pins changed (afc-data's own_to_host2979Fc.py), on A-fc's
# allocation folder (its FILL and PROFILE; ORDER, DISPATCH, PROVIDERS are the community's #2979 files, byte for byte)
for f in ORDER.json DISPATCH.json PROVIDERS.json; do cmp -s $AL/$f $D/allocation/$f || { echo "A-fc's $f is not #2979's"; bad=1; }; done
[ "$(shaf $AL/FILL.json.gz)" = 7fc766c6c3e1561cee60b258138a918020b96d85ea782a2b01ad0ee035fd082e ] && echo "   A-fc's FILL 7fc766c6..." || { echo "A-FC FILL DIFFERS"; bad=1; }
[ "$(shaf $AL/PROFILE.json)" = 5df7941120150492c64e62dc3362095eff3dfac9facff041ebb5c25b2e78f75e ] && echo "   A-fc's PROFILE 5df79411..." || { echo "A-FC PROFILE DIFFERS"; bad=1; }
cmp -s $AL/FILL.json.gz $PKG/fill/FILL.json.gz && cmp -s $AL/PROFILE.json $PKG/fill/PROFILE.json && echo "   = ../fill/FILL.json.gz and ../fill/PROFILE.json" || { echo "TREE FILL IS NOT ../fill/'S"; bad=1; }
mkdir -p $O/afd/out $O/afd/rerun/disp $O/afd/in; cp -R $AT $O/afd/tools; cp -R $AL $O/afd/in/alloc2979Fc
run $Q $O/afd/tools/own_to_host2979Fc.py $O/afd/in/alloc2979Fc $D/data $V $O/tables/residuals-2979.json $O/afd/out/H2979Fc
[ "$(shaf $O/afd/out/H2979Fc-hosttab.json.gz)" = a4baa4e69f2b7d0ff633fd176039d2759399c3c97a6176e8b0a9b2646b78963e ] \
  && echo "   host table H2979Fc-hosttab.json.gz: a4baa4e6... (the one the Lean data was made from)" || { echo "HOST TABLE DIFFERS"; bad=1; }
run $Q $O/afd/tools/fill_indep2979Fc.py $O/afd/out/H2979Fc-hosttab.json.gz $O/tables/H2979-hosttab.json.gz $O/afd/out/fill-independence.json
[ "$(shaf $O/afd/out/fill-independence.json)" = 955534d08ca47275a8343a4cd6b109bfd3749fa23d5e2ab41370304b342a3579 ] \
  && echo "   fill independence: every call field of the A-fc table = #2979's table (fill-independence.json 955534d0...)" || { echo "FILL INDEPENDENCE DIFFERS"; bad=1; }
# the head count (16,576 = 8,288 + 8,288 literal heads; no other rank-48 block): heads2979Fc.py prints ALL OK or SOME
# CHECK FAILED (exit 0 either way), so the verdict line is checked; one head more, or one other rank-48 block, must fail
$Q $O/afd/tools/heads2979Fc.py $O/afd/out/H2979Fc-hosttab.json.gz $AL/FILL.json.gz $V/TEMPLATE.json.gz - \
    $O/afd/out/heads2979Fc.json.gz $O/afd/out/heads2979Fc-summary.json 8288,8288,0,0 > $O/tables/heads.out 2>&1
tail -1 $O/tables/heads.out | grep -qx 'ALL OK' && echo "   head count 16,576 (8,288 + 8,288), no other rank-48 block: ALL OK" || { echo "HEAD COUNT CHECK FAILED"; bad=1; }
for n in 8289,8288,0,0 8288,8288,1,0; do
  $Q $O/afd/tools/heads2979Fc.py $O/afd/out/H2979Fc-hosttab.json.gz $AL/FILL.json.gz $V/TEMPLATE.json.gz - \
      $O/tables/never-heads-$n.json.gz $O/tables/never-heads-$n-summary.json $n > $O/tables/heads-neg-$n.out 2>&1
  tail -1 $O/tables/heads-neg-$n.out | grep -qx 'SOME CHECK FAILED' && echo "refused as required: counts $n" || { echo "COUNTS $n NOT REFUSED"; bad=1; }; done
run $Q $O/afd/tools/mk_parts.py $O/afd/out/H2979Fc-hosttab.json.gz H2979Fc 25 $O/lean/Work/XHostT/Data/H2979Fc
mkdir -p $O/neg/Work/XHostT/Data/H2979Fc $O/side/Work/XHostT/Data
mv $O/lean/Work/XHostT/Data/H2979Fc/NegGlob.lean $O/neg/Work/XHostT/Data/H2979Fc/NegGlob.lean   # a negative control: compared with ../controls/
refuse 3 $Q $O/afd/tools/mk_parts.py $O/afd/out/H2979Fc-hosttab.json.gz H2979Fcnload 25 $O/never-load --perturb load
refuse 3 $Q $O/afd/tools/mk_parts.py $O/afd/out/H2979Fc-hosttab.json.gz H2979Fcnwin 25 $O/never-win --perturb win
# the same perturbations with --emit-neg: the perturbed Calls modules (side files) and the NegSome / NegNone controls
for k in load win; do run $Q $O/afd/tools/mk_parts.py $O/afd/out/H2979Fc-hosttab.json.gz H2979Fcn$k 25 $O/neg/Work/XHostT/Data/H2979Fcn$k --perturb $k --emit-neg
  mkdir -p $O/side/Work/XHostT/Data/H2979Fcn$k; mv $O/neg/Work/XHostT/Data/H2979Fcn$k/Calls.lean $O/side/Work/XHostT/Data/H2979Fcn$k/Calls.lean; done
mkdir -p $O/tables/disp $O/lean/Work/XHostT/Data/D2979
run $Q $O/afd/tools/mk_disp2979.py $AL/DISPATCH.json $AL/ORDER.json 74 416 250 $O/tables/disp/Disp.lean
cp $O/tables/disp/Disp.lean $O/lean/Work/XHostT/Data/D2979/Disp.lean; cp $O/tables/disp/Disp.lean $O/afd/rerun/disp/Disp.lean
mkdir -p $O/never-disp; refuse 3 $Q $O/afd/tools/mk_disp2979.py $AL/DISPATCH.json $AL/ORDER.json 74 417 249 $O/never-disp/Disp.lean

echo "4. the instance modules (afc-data's checked renaming of the #2979 modules, numbers from the A-fc data above)"
run $Q $O/afd/tools/mk_inst2979Fc.py $O/lean
for f in NegHeads NegF1; do mv $O/lean/Work/XHostT/Data/H2979Fc/$f.lean $O/neg/Work/XHostT/Data/H2979Fc/$f.lean; done   # negative controls
mkdir -p $O/side/Work/XHostT/Inst
for f in S2979FcHB S2979FcHFW S2979FcUnit; do mv $O/lean/Work/XHostT/Inst/$f.lean $O/side/Work/XHostT/Inst/$f.lean; done   # side checks, not in the closure

echo "5. the FULL cut's profile, rate fact, statement chain, warp heads-first files, providers and the controls"
mkdir -p $O/afc
run $Q $CT/mk_profile2979Fc.py $AL/PROFILE.json $O/afc/profile2979Fc.json
cmp -s $O/afc/profile2979Fc.json $CT/profile2979Fc.json && echo "   profile2979Fc.json regenerated byte-identical (16,576 heads removed, target 25; nothing split; largest other rank 23)" || { echo "PROFILE2979FC DIFFERS"; bad=1; }
refuse 2 $Q $CT/mk_profile2979Fc.py $D/allocation/PROFILE.json $O/afc/never-profile.json
run $Q $CT/mk_rate_mixed.py scan 2979Fc=$O/afc/profile2979Fc.json
grep -E "17887045|17887046|largest|refus" $LOG | tail -4 | sed 's/^/   scan: /' | cut -c1-200
run $Q $CT/mk_rate_mixed.py emit 2979Fc=$O/afc/profile2979Fc.json 17887045 $O/lean
mkdir -p $O/never-rate; refuse 2 $Q $CT/mk_rate_mixed.py emit 2979Fc=$O/afc/profile2979Fc.json 17887046 $O/never-rate
mkdir -p $O/side/Work/Mixed; mv $O/lean/Work/Mixed/Rate2979FcSharp.lean $O/side/Work/Mixed/Rate2979FcSharp.lean   # a side check (../controls/side/)
run $Q $CT/mk_l2979Fc.py final $O/final-L2979Fc
refuse 2 $Q $CT/mk_l2979Fc.py final $O/never-L2979Fc --num 17887046
refuse 2 $Q $CT/mk_l2979Fc.py hyp $O/never-L2979Fch --num 17887046
cp -R $O/final-L2979Fc/stage/. $O/lean/
run $Q $X/x-inst/tools/mk_hfw2979.py $L/p-bridge/lean2 $O/lean
run $Q $CT/mk_prov2979Fc.py $O/lean $O/lean/Work/XHostT/Inst/S2979FcFacts.lean
run $Q $FT/mk_controls2979Fc.py $O/controls-made $O/side/Work/Mixed/Rate2979FcSharp.lean $O/neg
for f in Heads2979FcNeg.lean Rate2979FcSharpCtl.lean Ctl2979Fc.lean NegHeads.lean NegGlob.lean NegF1.lean \
         perturbed/NegSome119.lean perturbed/NegNone119.lean perturbed/NegSome46.lean perturbed/NegNone46.lean perturbed/NegSome334.lean perturbed/NegNone334.lean; do
  cmp -s $O/controls-made/$f $PKG/controls/$f && echo "   control $f regenerated byte-identical" || { echo "CONTROL $f DIFFERS"; bad=1; }; done
for f in Work/Mixed/Rate2979FcSharp.lean Work/XHostT/Inst/S2979FcHB.lean Work/XHostT/Inst/S2979FcHFW.lean Work/XHostT/Inst/S2979FcUnit.lean \
         Work/XHostT/Data/H2979Fcnload/Calls.lean Work/XHostT/Data/H2979Fcnwin/Calls.lean; do
  cmp -s $O/side/$f $PKG/controls/side/$f && echo "   side check $f regenerated byte-identical" || { echo "SIDE CHECK $f DIFFERS"; bad=1; }; done

echo "6. comparison with ../SHA256SUMS"
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
