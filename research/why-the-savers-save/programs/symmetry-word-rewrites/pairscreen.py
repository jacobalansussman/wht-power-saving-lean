"""Re-association screen: two reads into the same register t (t += a x ; t += b y, consecutive on t's chain, movable
together) replaced by one read t += c z of a register z that holds, at that moment, the clean row a x + b y (up to a
dyadic factor).  Frame-only screen with the local MILP; dirt repair considered only for winners."""
import sys, importlib.util, pickle, time
sys.path.insert(0,'.')
from fractions import Fraction as F
from collections import Counter, defaultdict
spec = importlib.util.spec_from_file_location('al', 'audit_local_py39.py'); al = importlib.util.module_from_spec(spec); spec.loader.exec_module(al)
import fw, search, reorder
w = fw.load('K8-FULL-WORD.json.gz'); v=112; R=480; n=2*v+R
W = search.Word(w); O = reorder.OW(W)
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
clean = lambda d: {k: x for k, x in d.items() if k < 2*v}
def key(d):
    if not d: return None, None
    k0 = min(d); c = d[k0]
    return tuple(sorted((k, x / c) for k, x in d.items())), c
nA = len(w['A']); nC = len(w['cut'])
gid_of = lambda idx: idx if idx < nA else idx - nC
# snapshot of clean rows before each main gate is too big; process online: for each gate index i, find partner j
# (next gate on t's chain) and test at time i.
nextread = {}
for r, ks in W.pos.items():
    for a, b in zip(ks, ks[1:]): nextread[a] = b
cands = []
hk = {}; hc = {}; groups = defaultdict(set)
for r in range(n):
    hk[r], hc[r] = key(clean(rows[r])); groups[hk[r]].add(r)
for idx, e in enumerate(main):
    if e[0] == 'add' and not (nA <= idx < nA + nC):
        g = gid_of(idx); t = e[1]
        h = nextread.get(g)
        if h is not None and O.sec[h] == O.sec[g]:
            e2 = W.seq[h][2]
            if e2[0] == 'add' and e2[1] == t and e2[2] != e[2] and O.commute(g, h):
                a = F(e[3], e[4]); b = F(e2[3], e2[4]); x, y = e[2], e2[2]
                tgt = defaultdict(F)
                for k, c in clean(rows[x]).items(): tgt[k] += a * c
                for k, c in clean(rows[y]).items(): tgt[k] += b * c
                tgt = {k: c for k, c in tgt.items() if c}
                kk, cc = key(tgt)
                if kk is not None:
                    for z in groups.get(kk, ()):
                        if z in (t, x, y): continue
                        q = cc / hc[z]
                        if q.denominator & (q.denominator - 1) == 0: cands.append((g, h, t, x, y, z, q))
    t = e[1]; old = hk[t]
    al.apply(e, rows)
    groups[old].discard(t); hk[t], hc[t] = key(clean(rows[t])); groups[hk[t]].add(t)
print('re-association options', len(cands), Counter('tgt' if v <= c[2] < 2*v else 'other' for c in cands), flush=True)
res = Counter(); good = []
for (g, h, t, x, y, z, q) in cands:
    # need h movable to right after g: all between commute with h
    pg, ph = O.posof[g], O.posof[h]
    if not all(O.commute(h, k) for k in O.order[pg+1:ph]): res['h not movable'] += 1; continue
    oldregs_g = O.regs[g]; oldorder = list(O.order); oldP, oldQ = O.P, O.Q
    O.regs[g] = (t, z); O.order.remove(h)
    O.refresh()
    d = None
    if O.feasible:
        changed = {k for k in O.order if O.P[k] != oldP.get(k) or O.Q[k] != oldQ.get(k)} | {g}
        r = O.local_milp(changed)
        if r is not None:
            newf, obj, regsF = r
            regsF = set(regsF) | {t, x, y, z}
            newfr = dict(O.f); newfr.update(newf)
            newc = sum(O.reg_cost(q_, newfr) for q_ in regsF)
            O.regs[g] = oldregs_g; O.order = oldorder; O.refresh_chains_only()
            oldc = sum(O.reg_cost(q_, O.f) for q_ in regsF)
            d = newc - oldc
    O.regs[g] = oldregs_g; O.order = oldorder; O.refresh()
    if d is None: res['infeasible'] += 1
    else:
        res[round(d, 3)] += 1
        if d < -1e-9: good.append((d, g, h, t, x, y, z, q))
print(sorted(res.items(), key=lambda x: str(x[0])))
good.sort(); print(good[:20])
pickle.dump(dict(cands=cands, good=good), open('pairs.pkl', 'wb'))
