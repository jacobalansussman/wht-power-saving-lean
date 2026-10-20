#!/bin/bash
# REGENERATE.sh OUTDIR: make every fill-dependent generated file of the design-B lever instance (L3945H: #3945 with
# every head kept) again, from the community's #3945 owner folder (tree/wht45/lever/g-heads/dl/pr3945/owner, PR head
# 91f81d13) with the design's fill laid over it (../fill/), with the programs in tree/ (the workspace's own files at
# the same relative paths; the only change is the path-only patch in PATCHES.txt), and compare the result byte for
# byte with ../SHA256SUMS (overlay/...) and with the controls in ../controls/.  Each generator that writes a figure, a
# count or a table is also run one step too far and must refuse.  Light Python 3 (standard library) only; about
# 4 minutes and under 1 GB.  bB-final (wht45/allheads, 19 Oct 2026), x-final's REGENERATE.sh with the design-B steps.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); PKG=$(dirname "$HERE")
export GEN_ROOT=$HERE/tree
O=${1:?usage: REGENERATE.sh OUTDIR}; [ -e "$O" ] && { echo "REFUSE: $O exists"; exit 1; }; mkdir -p "$O"; O=$(cd "$O" && pwd)
Q="python3 -I -B"
L=$GEN_ROOT/wht45/lever; BC=$GEN_ROOT/wht43/bound-check; BT=$GEN_ROOT/wht45/allheads/bB-data/tools
GH=$L/g-heads; V=$GH/dl/pr3945/owner
LOG=$O/regenerate.log; : > $LOG
bad=0
run() { echo "+ $*" >> $LOG; "$@" >> $LOG 2>&1 || { echo "FAILED: $*"; bad=1; }; }
refuse() { local want=$1; shift; echo "+ (must exit $want) $*" >> $LOG; "$@" >> $LOG 2>&1; local rc=$?
  if [ $rc = $want ]; then echo "refused as required (exit $rc): ${*: -4}"; else echo "DID NOT REFUSE (exit $rc): $*"; bad=1; fi; }
shaf() { shasum -a 256 "$1" | cut -c1-64; }
(cd $HERE/tree && shasum -a 256 -c --quiet ../tree.sha256) || { echo "generators tree differs from tree.sha256"; exit 1; }

echo "0. the design's owner folder: the community's #3945 owner folder with the design-B fill laid over it"
cp -R $V $O/owner
cp $PKG/fill/FILL.json.gz $O/owner/FILL.json.gz
cp $PKG/fill/PROFILE.json $O/owner/PROFILE.json
cp $PKG/fill/old-FILL.json.gz $O/owner/inputs/old-FILL.json.gz
(cd $O/owner && shasum -a 256 -c --quiet $HERE/owner-design.sha256) && echo "   owner folder = the one the Lean data was made from (24 files; FILL 9c06677d, old-FILL 0e61fce2)" \
  || { echo "OWNER FOLDER DIFFERS"; bad=1; }

echo "1. the head gate G0 in our lowering (g-heads' g0_staged.py)"
run $Q $BT/g0_staged.py $BC $O/owner $O/owner/FILL.json.gz $GH/data/gcert/S3945-base.gcert1b.json.gz \
    $GH/data/gcert/S3945-extra.gcert1b.json.gz $O/owner/PROFILE.json $O/g0-3945H.json --label "#3945H (all heads)"
$Q - $O/g0-3945H.json <<'PY' || bad=1
import json, sys
d = json.load(open(sys.argv[1])); g = d['G0']
ok = (d['ok'] and g['rank48_by_class'] == {'HEAD-literal': 80640, 'neither': 57312} and g['violations'] == {}
      and g['cross_replica_roles'] == 1161423 and d['call_blocks']['global_profile'] == {'1': 8041680, '2': 2306880, '3': 729648, '4': 353160, '5': 69120, '8': 25920, '48': 137952})
print('   G0: rank-48 bins', g['rank48_by_class'], '; window violations', len(g['violations']), '; profile = #3945\'s:', ok)
sys.exit(0 if ok else 1)
PY

echo "2. the host table of the fill (p-own's converter), its head list and its Lean parts (p-hostck's generator)"
mkdir -p $O/tables $O/lean/Work/XHostT/Data
run $Q $BT/own_to_host3945H.py $O/owner $GH/data/residuals-3945.json $O/g0-3945H.json $O/tables/S5H
[ "$(shaf $O/tables/S5H-hosttab.json.gz)" = b3185813426853705e07c8d374b116252e4e4576e7f1f8dcf6069e8d239dfbf1 ] \
  && echo "   host table S5H-hosttab.json.gz: b3185813... (the one the Lean data was made from)" || { echo "HOST TABLE DIFFERS"; bad=1; }
$Q $BT/heads3945H.py $O/tables/S5H-hosttab.json.gz $O/owner/FILL.json.gz $O/owner/inputs/TEMPLATE.json.gz $O/g0-3945H.json \
    $O/tables/heads3945H.json.gz $O/tables/heads3945H-summary.json > $O/tables/heads.out 2>&1
tail -1 $O/tables/heads.out | grep -qx 'ALL OK' && echo "   head list: 80,640 literal heads (40,320 + 40,320), 57,312 other rank-48 bins: ALL OK" || { echo "HEAD LIST CHECK FAILED"; bad=1; }
run $Q $BT/mk_parts.py $O/tables/S5H-hosttab.json.gz S5H 25 $O/lean/Work/XHostT/Data/S5H
mkdir -p $O/side; mv $O/lean/Work/XHostT/Data/S5H/NegGlob.lean $O/side/NegGlob.lean   # a negative control (must fail): ../controls/
mv $O/lean/Work/XHostT/Data/S5H/MANIFEST.json $O/side/S5H-MANIFEST.json
refuse 3 $Q $BT/mk_parts.py $O/tables/S5H-hosttab.json.gz S5Hn 25 $O/never-load --perturb load
refuse 3 $Q $BT/mk_parts.py $O/tables/S5H-hosttab.json.gz S5Hw 25 $O/never-win --perturb win

echo "3. the dispatch (p-disp's generator): it does not depend on the fill, so Disp.lean must be #3945's byte for byte"
mkdir -p $O/tables/disp $O/lean/Work/XHostT/Data/D3945
run $Q $L/p-disp/tools/mk_disp3945.py $O/owner/inputs/DISPATCH.json $O/owner/inputs/ORDER.json 360 1672 1568 $O/tables/disp/Disp.lean
cp $O/tables/disp/Disp.lean $O/lean/Work/XHostT/Data/D3945/Disp.lean      # (its side file D3945-disp.json stays in tables/disp/)
mkdir -p $O/never-disp; refuse 3 $Q $L/p-disp/tools/mk_disp3945.py $O/owner/inputs/DISPATCH.json $O/owner/inputs/ORDER.json 360 1673 1567 $O/never-disp/Disp.lean

echo "4. the rate fact at the figure (e-r45m1's generator) and its sharpness control"
run $Q $BT/mk_rate_mixed.py scan 3945H=$BT/profile3945H.json
run $Q $BT/mk_rate_mixed.py emit 3945H=$BT/profile3945H.json 17581520 $O/lean
mkdir -p $O/never-rate/lean; refuse 2 $Q $BT/mk_rate_mixed.py emit 3945H=$BT/profile3945H.json 17581521 $O/never-rate/lean
mv $O/lean/Work/Mixed/Rate3945HSharp.lean $O/side/Rate3945HSharp.lean    # in ../controls/Rate3945HSharpCtl.lean

echo "5. the hand-written fill-dependent files, instantiated for the fill (bB-data's mk_3945H.py: exact counted replacements)"
run $Q $BT/mk_3945H.py $O/inst --neg-heads
for f in NegHeads NegHeadsOld NegRes; do mv $O/inst/Work/XHostT/Data/S5Hneg/$f.lean $O/side/$f.lean; done
rmdir $O/inst/Work/XHostT/Data/S5Hneg
mv $O/inst/Work/XHostT/Inst/S5HHB.lean $O/side/S5HHB.lean           # a side check (80,640 head bins), not in the closure
cp -R $O/inst/. $O/lean/

echo "6. the statement chain at the figure (i-skel's mk_l3945.py with the design's profile)"
run $Q $BT/mk_l3945.py final $O/final-L3945H --suffix L3945H --num 17581520 --profile 3945H=$BT/profile3945H.json
refuse 2 $Q $BT/mk_l3945.py final $O/never-L3945H --suffix L3945H --num 17581521 --profile 3945H=$BT/profile3945H.json
cp -R $O/final-L3945H/stage/. $O/lean/

echo "7. comparison with ../SHA256SUMS (overlay) and ../controls/"
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
    for p in extra: f.write('made, not in the overlay\t' + p + '\n')
    for p in notmade: f.write('in the overlay, not generated here (fill-independent, as in the #3945 build: see README.txt)\t' + p + '\n')
print(f'   regenerated overlay files byte-identical: {len(same)}; DIFFERENT: {len(diff)}; made but not in the overlay: {len(extra)}; '
      f'overlay files not generated here: {len(notmade)} (list: {out})')
for p in diff: print('DIFFERS', p)
for p in extra: print('made, not in the overlay:', p)
sys.exit(1 if diff or extra or not same else 0)
PY
for f in NegHeads NegHeadsOld NegRes NegGlob; do
  cmp -s $O/side/$f.lean $PKG/controls/$f.lean && echo "   controls/$f.lean: byte-identical" || { echo "CONTROL $f DIFFERS"; bad=1; }; done
$Q - $O/side/Rate3945HSharp.lean $PKG/controls/Rate3945HSharpCtl.lean <<'PY' || bad=1
import sys
gen, ctl = (open(f).read() for f in sys.argv[1:3])
body = '\n'.join(l for l in gen.split('\n') if not l.startswith('import '))
ok = ctl.endswith(body.lstrip('\n')) or body.strip() in ctl
print('   controls/Rate3945HSharpCtl.lean contains the regenerated Rate3945HSharp.lean below its imports:', ok)
sys.exit(0 if ok else 1)
PY
if [ $bad = 0 ]; then echo "REGENERATE: ALL AS EXPECTED (log $LOG)"; exit 0; else echo "REGENERATE: SOMETHING FAILED (log $LOG)"; exit 1; fi
