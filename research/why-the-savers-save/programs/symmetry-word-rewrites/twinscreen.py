"""Screen 'read the twin' rewrites: gate t += q*s at a moment when another register s2 holds the same clean row
(source/target columns; helper dirt may differ) up to a dyadic factor.  Frame-only screen: swap s -> s2 in the gate,
re-optimise frames by local MILP, report the cost delta.  Scalar repair (dirt) is checked later for any winner."""
import sys, importlib.util, pickle, time
sys.path.insert(0,'.')
from fractions import Fraction as F
from collections import Counter, defaultdict
spec = importlib.util.spec_from_file_location('al', 'audit_local_py39.py'); al = importlib.util.module_from_spec(spec); spec.loader.exec_module(al)
import fw, search, reorder
w = fw.load('K8-FULL-WORD.json.gz'); v=112; R=480; n=2*v+R
ident = lambda: [{k: F(1)} for k in range(n)]
main = w['A'] + w['cut'] + w['B']
raw = ident()
for e in main: al.apply(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2*v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2*v}
for e in w['precondition']: al.apply(e, dirt)
comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]
rows = ident()
for e in w['precondition'] + comp: al.apply(e, rows)
def key(d):
    d = {k: x for k, x in d.items() if k < 2*v}
    if not d: return None, None
    k0 = min(d); c = d[k0]
    return tuple(sorted((k, x / c) for k, x in d.items())), c
hk = {}; hc = {}
groups = defaultdict(set)
for r in range(n):
    hk[r], hc[r] = key(rows[r]); groups[hk[r]].add(r)
nA = len(w['A']); nC = len(w['cut'])
twins = []
for idx, e in enumerate(main):
    if e[0] == 'add' and not (nA <= idx < nA + nC) and hk[e[2]] is not None:
        s = e[2]
        for s2 in groups[hk[s]] - {s, e[1]}:
            fac = hc[s] / hc[s2]           # row(s) = fac * row(s2) on clean columns
            q = F(e[3], e[4]) * fac
            if q.denominator & (q.denominator - 1) == 0:
                gid = idx if idx < nA else idx - nC
                twins.append((gid, s, s2, q))
    t = e[1]; old = hk[t]
    al.apply(e, rows)
    groups[old].discard(t); hk[t], hc[t] = key(rows[t]); groups[hk[t]].add(t)
print('dyadic twin options', len(twins), flush=True)
W = search.Word(w); O = reorder.OW(W)
res = Counter(); good = []
t0 = time.time()
for (g, s, s2, q) in twins:
    oldregs = O.regs[g]
    newregs = (oldregs[0], s2)
    oldP, oldQ = O.P, O.Q
    O.regs[g] = newregs; O.W.regs = O.regs
    O.refresh()
    ok = O.feasible
    d = None
    if ok:
        changed = {h for h in O.order if O.P[h] != oldP.get(h) or O.Q[h] != oldQ.get(h)} | {g}
        r = O.local_milp(changed)
        if r is not None:
            newf, obj, regsF = r
            regsF = set(regsF) | {s, s2}
            newfr = dict(O.f); newfr.update(newf)
            newc = sum(O.reg_cost(x, newfr) for x in regsF)
            O.regs[g] = oldregs; O.refresh_chains_only()
            oldc = sum(O.reg_cost(x, O.f) for x in regsF)
            d = newc - oldc
    O.regs[g] = oldregs; O.refresh()
    if d is None: res['infeasible'] += 1
    else:
        res[round(d, 3)] += 1
        if d < -1e-9: good.append((d, g, s, s2, q))
print(sorted(res.items(), key=lambda x: str(x[0])), 'secs', round(time.time()-t0,1))
good.sort()
for x in good[:30]: print(x)
pickle.dump(dict(twins=twins, good=good), open('twins.pkl', 'wb'))
