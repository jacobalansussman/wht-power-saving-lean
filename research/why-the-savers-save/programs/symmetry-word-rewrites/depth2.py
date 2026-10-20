"""Depth-2 search over NON-cosmetic commuting reorders: every single chain-adjacent reorder that changes some frame
interval (P or Q) is evaluated (local MILP); then, for each with d1 < D1MAX, every non-cosmetic reorder on the registers
it disturbed is tried on top.  Reports the best d1+d2."""
import sys, os, time, json
sys.path.insert(0, '.')
import fw, search, reorder
from collections import Counter
D1MAX = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz'))); O = reorder.OW(W)
def cand_moves(regset=None):
    mv = []
    for r, gs in O.chain.items():
        if regset is not None and r not in regset: continue
        for g1, g2 in zip(gs, gs[1:]):
            if O.sec[g1] != O.sec[g2] or not O.commute(g1, g2): continue
            i, j = O.posof[g1], O.posof[g2]; between = O.order[i + 1:j]
            if all(O.commute(g2, h) for h in between): mv.append((g2, i))
            if all(O.commute(g1, h) for h in between): mv.append((g1, j))
    return list(dict.fromkeys(mv))
def evaluate(g, p):
    """apply move tentatively; return (d, touched registers) or None; always undoes."""
    oldorder = list(O.order); oldP, oldQ = O.P, O.Q
    O.order.remove(g); O.order.insert(p, g); O.refresh()
    res = None
    if O.feasible:
        changed = {h for h in O.order if O.P[h] != oldP[h] or O.Q[h] != oldQ[h]}
        if changed:
            r = O.local_milp(changed | {g})
            if r is not None:
                newf, obj, regsF = r
                newfr = dict(O.f); newfr.update(newf)
                newc = sum(O.reg_cost(q, newfr) for q in regsF)
                O.order = oldorder; O.refresh_chains_only()
                oldc = sum(O.reg_cost(q, O.f) for q in regsF)
                touched = {q for h in changed for q in O.regs[h]} | set(O.regs[g])
                res = (newc - oldc, touched, newf)
    O.order = oldorder; O.refresh()
    return res
t0 = time.time()
base = cand_moves(); first = []
for g, p in base:
    r = evaluate(g, p)
    if r is not None: first.append((r[0], g, p, r[1]))
print('single moves', len(base), 'non-cosmetic feasible', len(first), 'best d1', min(x[0] for x in first), Counter(round(x[0], 2) for x in first).most_common(8), f'{time.time()-t0:.0f}s', flush=True)
best = (0.0, None); n2 = 0
for d1, g, p, touched in sorted(first, key=lambda x: x[0]):
    if d1 >= D1MAX: continue
    snap_order = list(O.order); snap_f = dict(O.f)
    # apply m1 for real (with its re-framing)
    r = O.try_move(g, p, accept=lambda d: True)
    if r is None: O.order = snap_order; O.f = snap_f; O.refresh(); continue
    for g2, p2 in cand_moves(touched):
        if (g2, p2) == (g, p): continue
        r2 = evaluate(g2, p2)
        if r2 is None: continue
        n2 += 1
        tot = r + r2[0]
        if tot < best[0] - 1e-9:
            best = (tot, ((g, p), (g2, p2)))
            print('IMPROVING pair', round(tot, 4), best[1], flush=True)
    O.order = snap_order; O.f = snap_f; O.refresh()
print('pairs evaluated', n2, 'best', best, f'{time.time()-t0:.0f}s', flush=True)
json.dump(dict(single=len(base), noncosmetic=len(first), pairs=n2, best=best[0], moves=best[1]), open('depth2_%s.json' % os.path.basename(os.environ.get('WORD', 'k8')).split('.')[0], 'w'))
