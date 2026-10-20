# c3_s8.py -- make the singular subgame explicit as S8: points, h, K, B', the affine structure; check Young reps.
import numpy as np, itertools
from collections import Counter
from lib514 import *
from libs8 import *
print('each nonsingular vector in exactly two point-cliques:', Counter(len(pair_of(v)) for v in NS))
refl = [TRANS[v] for v in NS]
Gp = gen_closure(refl)
perms = {tuple(to_perm(g)) for g in Gp}
print('|O(q0)| =', len(Gp), ' distinct images in S8:', len(perms))
rng = np.random.default_rng(0)
ok = True
for _ in range(200):
    a = Gp[rng.integers(len(Gp))]; b = Gp[rng.integers(len(Gp))]
    ok &= np.array_equal(to_perm(comp(a, b)), compose(to_perm(a), to_perm(b)))
print('to_perm is a homomorphism (200 random pairs):', ok)
ph = to_perm(h1); pK = to_perm(K)
def cyc(p):
    seen = set(); out = []
    for i in range(8):
        if i in seen: continue
        c = [i]; seen.add(i); j = int(p[i])
        while j != i: c.append(j); seen.add(j); j = int(p[j])
        if len(c) > 1: out.append(tuple(c))
    return out
print('h ->', cyc(ph), ' h2 ->', cyc(to_perm(h2)), ' h3 ->', cyc(to_perm(h3)), ' K ->', cyc(pK))
Bp = [g for g in Gp if frozenset(int(g[v]) for v in LZ) == LZ]
PB = [to_perm(g) for g in Bp]
print('|B\'| =', len(PB))
# B' orbits on k-subsets
for k in (1, 2, 3, 4):
    subs = list(itertools.combinations(range(8), k)); seen = set(); orbs = []
    for s in subs:
        if s in seen: continue
        o = {tuple(sorted(int(p[i]) for i in s)) for p in PB}; seen |= o; orbs.append(len(o))
    print(f"B' orbits on {k}-subsets:", sorted(orbs))
blocks = None
subs = list(itertools.combinations(range(8), 4)); seen = set()
for s in subs:
    if s in seen: continue
    o = {tuple(sorted(int(p[i]) for i in s)) for p in PB}; seen |= o
    if len(o) == 14: blocks = sorted(o)
print('the 14 affine planes (blocks of the Steiner system S(3,4,8) of L_Z):', blocks)
# Levi image
PL = [to_perm(l) for l in LEVI]
print('Levi fixes points:', sorted({i for i in range(8) if all(p[i] == i for p in PL)}),
      ' Levi orbits on points:', sorted(Counter(min(int(p[i]) for p in PL) for i in range(8)).values()))
# K's affine structure
blocksK = sorted({tuple(sorted(int(pK[i]) for i in b)) for b in blocks})
print('blocks of K(L_Z) structure:', blocksK, ' shared blocks with L_Z:', len(set(blocks) & set(blocksK)))
# Young reps: homomorphism and dimension checks
dims = {}
for lam in partitions(8):
    R = YoungRep(lam); dims[lam] = R.dim
    err = 0
    for _ in range(20):
        a = Gp[rng.integers(len(Gp))]; b = Gp[rng.integers(len(Gp))]
        pa, pb = to_perm(a), to_perm(b)
        err = max(err, np.abs(R(compose(pa, pb)) - R(pa) @ R(pb)).max())
    tr_t = np.trace(R(ph))
    print(lam, 'dim', R.dim, ' hom err %.1e' % err, ' trace(transposition)=%.0f' % tr_t, ' trace(K)=%.0f' % np.trace(R(pK)))
print('sum of squares of dims:', sum(d * d for d in dims.values()))
np.save('s8_points_h_K.npy', np.array([ph, pK]))
