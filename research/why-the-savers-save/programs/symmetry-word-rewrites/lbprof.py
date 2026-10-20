"""Profiles (climb-size counts) of the per-register decoupled optima, for pricing the two ceilings as what-ifs:
 (a) fixed order, intervals [P,Q] (valid ceiling for any frame-only change of this order);
 (b) dependency-only intervals [Pd,Qd] (optimistic gauge for reorders; not a valid ceiling)."""
import sys, json; sys.path.insert(0,'.')
import fw, search, ilp, reorder2
from search import le, FC
from ilp_wide import all_between
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz')); O = reorder2.OW2(W)
P, Q = ilp.bounds_PQ(W)
d = json.load(open('lbdep.json')); Pd = {int(k): tuple(v) for k, v in d['Pd'].items()}; Qd = {int(k): tuple(v) for k, v in d['Qd'].items()}
def dp(r, lo, hi, chain):
    states = {W.start[r]: (0.0, Counter())}
    for k in chain:
        gap = len(hi[k]) - len(lo[k])
        cands = all_between(lo[k], hi[k]) if gap <= 3 else {lo[k], hi[k], W.f[k]}
        new = {}
        for U in cands:
            best = None
            for V, (c, p) in states.items():
                if le(V, U):
                    cc = c + (FC[len(U)-len(V)] if U != V else 0.0)
                    if best is None or cc < best[0]:
                        pp = Counter(p)
                        if U != V: pp[len(U)-len(V)] += 1
                        best = (cc, pp)
            if best is not None: new[U] = best
        states = new
    E = W.end[r]; best = None
    for V, (c, p) in states.items():
        if le(V, E):
            cc = c + (FC[len(E)-len(V)] if V != E else 0.0)
            if best is None or cc < best[0]:
                pp = Counter(p)
                if V != E: pp[len(E)-len(V)] += 1
                best = (cc, pp)
    return best
out = {}
for name, lo, hi, chainf in (('fixed_order_PQ', P, Q, lambda r: W.pos.get(r, [])), ('dependency_only', Pd, Qd, lambda r: O.chain.get(r, []))):
    tot = 0.0; prof = Counter()
    for r in range(W.n):
        c, p = dp(r, lo, hi, chainf(r)); tot += c; prof.update(p)
    prof[8] += 8   # centre paths unchanged
    out[name] = dict(cost=tot, profile={str(k): v for k, v in sorted(prof.items())})
    print(name, round(tot, 3), dict(sorted(prof.items())))
json.dump(out, open('lbprof.json', 'w'), indent=1)
