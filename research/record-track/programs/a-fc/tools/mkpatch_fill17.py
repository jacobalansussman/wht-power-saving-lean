# fc-reopt: fill17_fc.py = #2979's fill17.py + literal-head variables + the separated-cut objective (our derivative; the
# hub file is untouched). Head part = o-mip's tools/mkheads.py (itself o-model's probe/mkpatch.py), unchanged in effect:
#   --force-heads    : X[g(o, 8/10), empty] == 224 for every owner (H = 224 T)
#   --max-heads      : maximize H at the fixed global mix
#   --head-bonus K   : with --free-mix, subtract K (objective units) per head
#   --bstar B        : with --free-mix, linearise at B instead of the relaxation cap root
# New here:
#   --fc-C C         : with --free-mix, the bin cost is 73^b phi_C(r) (phi_C: every paid block above C split into
#                      floor(r/C) blocks of C and one of r mod C) instead of H(r, b) = 73^b r^(1-b): the moment of the
#                      separated cut of cap C (C = 24: the full cut), in the hub's H units
#   --catalogue file:PATH : with --free-mix, keep exactly the patterns listed in PATH ([[cap, counts], ...])
import sys
src = open(sys.argv[1]).read()
a = "    ap.add_argument('--method', default="
b = ("    ap.add_argument('--force-heads', action='store_true')\n    ap.add_argument('--max-heads', action='store_true')\n"
     "    ap.add_argument('--head-bonus', type=int, default=0)\n    ap.add_argument('--bstar', type=float, default=None)\n"
     "    ap.add_argument('--fc-C', type=int, default=None)\n" + a)
assert src.count(a) == 1; src = src.replace(a, b)
a = "        hinted = {}\n"
b = ("        head_groups = [gi for gi, g in enumerate(G) if groups[g][0][1] in (8, 10)]\n"
     "        assert len(head_groups) == T and all(G[gi][0] == 48 and all(bb[1] in (8, 10) for bb in groups[G[gi]]) for gi in head_groups)\n"
     "        q0 = tuple(0 for _ in WID)\n"
     "        hv = [vs[Xv[gi, q0]] for gi in head_groups if (gi, q0) in Xv]\n"
     "        log('HEAD GROUPS', len(head_groups), 'with an empty-pattern variable', len(hv))\n"
     "        if a.force_heads:\n"
     "            for gi in head_groups:\n"
     "                model.Add(vs[Xv[gi, q0]] == len(groups[G[gi]]))\n"
     "        if a.max_heads:\n"
     "            model.Maximize(cp_model.LinearExpr.Sum(hv))\n"
     "        if a.free_mix and a.head_bonus:\n"
     "            model.Minimize(cp_model.LinearExpr.WeightedSum([v for v, c in obj] + hv, [c for v, c in obj] + [-a.head_bonus] * len(hv)))\n"
     "            log('HEAD BONUS', a.head_bonus, 'per head on', len(hv), 'variables')\n"
     "        hinted = {}\n")
assert src.count(a) == 1; src = src.replace(a, b)
a = "            bstar = float((cap_lo + cap_hi) / 2)\n"
b = ("            bstar = float((cap_lo + cap_hi) / 2) if a.bstar is None else a.bstar\n            log('FREE-MIX bstar', bstar, 'fc-C', a.fc_C)\n"
     "            def _fcH(r, b, C=a.fc_C):\n"
     "                if C is None: return X.Hfloat(r, b)\n"
     "                q, m = divmod(r, C)\n"
     "                return (73.0 ** b) * (q * C ** (1 - b) + (m ** (1 - b) if m else 0.0))\n")
assert src.count(a) == 1; src = src.replace(a, b)
a = "                c = 0 if r == 0 else int(round(a.objective_scale * X.Hfloat(r, bstar)))\n"
b = "                c = 0 if r == 0 else int(round(a.objective_scale * _fcH(r, bstar)))\n"
assert src.count(a) == 1; src = src.replace(a, b)
a = "        elif a.catalogue == 'all':\n"
b = ("        elif a.catalogue.startswith('file:'):\n"
     "            keep = {(int(cap), tuple(int(k) for k in counts)) for cap, counts in json.load(open(a.catalogue[5:]))}\n"
     "        elif a.catalogue == 'all':\n")
assert src.count(a) == 1; src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
print('patched', sys.argv[2])
