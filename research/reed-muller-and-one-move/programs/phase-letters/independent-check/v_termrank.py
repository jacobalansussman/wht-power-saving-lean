"""Own exact test (valid, time-order free): rank E[D_F, D_F] <= term rank of the positions not forced to zero by tight
supports (diagonal at live moves; pairs sharing a tight support). Exclusion if term rank < |D_F| - r_F (rm-theory 8.5).
Run on rm2-phase's 411 S+D+C+Z classes of 4 in 14 (dirs 1..14); compare with its tiling verdicts."""
import sys, re, numpy as np
sys.path.insert(0, '../ends-only-block4')
from tile_alt import target_matrices
from vb import *
mats, nonalt = target_matrices(4); m = 4; Qz = Qz_matrix(m); nf, nd = Qz.shape
def matching(rows, cols, ok):
    match = {}
    def aug(u, seen):
        for v in cols:
            if ok(u, v) and v not in seen:
                seen.add(v)
                if v not in match or aug(match[v], seen):
                    match[v] = u; return True
        return False
    return sum(aug(u, set()) for u in rows)
def test(dirs, tg):
    inst = [z - 1 for z in dirs]; n = len(inst)
    cnt = np.zeros(nd, dtype=np.int64)
    for d in inst: cnt[d] += 1
    Dq = cnt @ Qz.T
    r = pair_matrix(m, [np.array(mats[t]) for t in tg]).sum(1)
    assert (Dq >= r).all()
    sets = [frozenset(k for k, d in enumerate(inst) if Qz[f, d]) for f in range(nf)]
    tight = [sets[f] for f in range(nf) if Dq[f] == r[f] and r[f] >= 1]
    live = set().union(*tight) if tight else set()
    zero = set()
    for s in tight:
        for u in s:
            for v in s: zero.add((u, v))
    for f in range(nf):
        D = sorted(sets[f]); need = len(D) - r[f]
        if need <= 0: continue
        tr = matching(D, D, lambda u, v: (u, v) not in zero)
        if tr < need: return ('excluded', f, len(D), int(r[f]), tr)
    return ('passes',)
out = {}
for line in open('../ends-only-block4/surv_fe_all_414.txt'):
    a, b = line.strip().split('|'); d = [int(x) for x in a.split(',')]; t = [int(x) for x in b.split(',')]
    v = test(d, t); key = (v[0], sum(x >= nonalt for x in t) > 0); out[key] = out.get(key, 0) + 1
print("term-rank test on the 411 classes of 4 in 14 (verdict, has alternating target): count:", out)
