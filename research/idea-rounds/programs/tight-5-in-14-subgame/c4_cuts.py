# c4_cuts.py -- cut certificates.
# (a) calibration: V_56 under the full free group P0 (10752 elements): best cut must be 14/5 = 2.8 (the hook).
# (b) the singular subgame: free group B' = AGL(3,2) (1344 elements); every irreducible of S8 (Young), and V_56
#     restricted to O(q0): best cut ratio f(K)/f(h). If anything beats 2.8, no 5-in-14 lives in the subgame.
import numpy as np, time
from lib514 import *
from libs8 import *
from libcut import *
t0 = time.time(); rng = np.random.default_rng(7)
R = Irreps()
def closure_with_mats(gens, rep):
    el = {key(ID): (ID, np.eye(rep(ID).shape[0]))}; fr = [key(ID)]
    gm = [(g, rep(g)) for g in gens]
    while fr:
        nx = []
        for k in fr:
            x, M = el[k]
            for g, G in gm:
                y = comp(g, x); ky = key(y)
                if ky not in el: el[ky] = (y, G @ M); nx.append(ky)
        fr = nx
    return [v[1] for v in el.values()], [v[0] for v in el.values()]
# (a)
P0_gens = A_GENS + [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]
mats, els = closure_with_mats(P0_gens, R.r56)
print('|P0| =', len(mats), ' max matrix error vs direct', max(np.abs(mats[i] - R.r56(els[i])).max() for i in range(0, len(els), 997)))
comps = isotypic(mats, rng)
print('V_56 under P0: components (d, m):', sorted((c['d'], c['m']) for c in comps))
top, allr = best_cuts(comps, R.r56(h1), R.r56(K), extra=[R.r56(comp(h1, h2))])
for t in top: print('   ratio %.4f dimY %d cap %d demand %d (2-step %d) subset %s' % (t[0], t[1], t[2], t[3], t[5], t[4]))
print('[%.0fs]' % (time.time() - t0))
# (b) subgame free group B'
refl = [TRANS[v] for v in NS]
Gp = gen_closure(refl)
Bp = [g for g in Gp if frozenset(int(g[v]) for v in LZ) == LZ]
Bp_gens = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
           TRANS[Z(1) | Z(2)], TRANS[Z(1) | Z(3)], TRANS[Z(2) | Z(3)]]
# check: the three transvections t_(0|w), w of weight 2 -> alternating? t_(0|w): z += (x.w) w, S = w w^T has diagonal w -> NOT alternating.
# use instead the alternating radical elements found by brute force
alt = [g for g in gen_closure(A_GENS) if all(q0(int(g[v])) == q0(v) for v in range(NV))]
Bp_gens = Bp_gens[:2] + [a for a in alt if not np.array_equal(a, ID)][:3]
assert len(gen_closure(Bp_gens)) == 1344
mats_b, els_b = closure_with_mats(Bp_gens, R.r56)
comps_b = isotypic(mats_b, rng)
print("V_56 under B' (singular subgame): components (d, m):", sorted((c['d'], c['m']) for c in comps_b))
top, allr = best_cuts(comps_b, R.r56(h1), R.r56(K), extra=[R.r56(comp(h1, h2))], top=8)
for t in top: print('   ratio %.4f dimY %d cap %d demand %d (2-step %d) subset %s' % (t[0], t[1], t[2], t[3], t[5], t[4]))
print('[%.0fs]' % (time.time() - t0))
ph = to_perm(h1); pK = to_perm(K); ph2 = to_perm(comp(h1, h2))
PBg = [to_perm(g) for g in Bp_gens]
best_overall = []
for lam in partitions(8):
    Y = YoungRep(lam)
    if Y.dim == 1: continue
    # B' elements as matrices via closure on permutations
    el = {tuple(range(8)): np.eye(Y.dim)}; fr = [tuple(range(8))]
    gm = [(p, Y(p)) for p in PBg]
    while fr:
        nx = []
        for k in fr:
            M = el[k]
            for p, G in gm:
                q = tuple(int(x) for x in compose(p, np.array(k)))
                if q not in el: el[q] = G @ M; nx.append(q)
        fr = nx
    assert len(el) == 1344
    comps = isotypic(list(el.values()), rng)
    top, allr = best_cuts(comps, Y(ph), Y(pK), extra=[Y(ph2)], top=3)
    t = top[0] if top else None
    print(lam, 'dim', Y.dim, " B'-components (d,m):", sorted((c['d'], c['m']) for c in comps),
          ' best cut: ratio %.4f cap %d demand %d 2-step %d dimY %d' % (t[0], t[2], t[3], t[5], t[1]) if t else 'none')
    best_overall.append((t[0] if t else 0, lam))
best_overall.sort(reverse=True)
print('best cut ratios over S8 irreducibles:', [(round(r, 4), l) for r, l in best_overall[:6]])
print('[%.0fs]' % (time.time() - t0))
