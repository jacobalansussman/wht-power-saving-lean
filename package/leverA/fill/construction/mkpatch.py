# o-model probe: make fill17_heads.py = #2979's fill17.py + head variables (our derivative, hub file untouched).
#   --force-heads : every bin of the child-8/10 group of every owner is empty (X[g, empty] == 224)
#   --max-heads   : CP-SAT maximizes sum_o X[g(o, 8/10), empty] (the literal head count) at the fixed global mix
# The literal head count of a T37 fill is exactly that sum (MODEL.md, section 4); doubling to T74 doubles it.
import sys
src = open(sys.argv[1]).read()
a = "    ap.add_argument('--method', default="
b = "    ap.add_argument('--force-heads', action='store_true')\n    ap.add_argument('--max-heads', action='store_true')\n" + a
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
     "        hinted = {}\n")
assert src.count(a) == 1; src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
