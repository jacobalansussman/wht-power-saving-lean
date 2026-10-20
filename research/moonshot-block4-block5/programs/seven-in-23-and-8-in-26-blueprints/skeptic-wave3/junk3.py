"""raw-matrix confirmation of the d = 2 failure of the Theorem R' condition: t = l2 D(z2) l1 D(z1) a (an O2 output term,
shared move in O O F O F).  q_g = t^-1 c_g t; S = {g : K q_g K^-1 monomial}; junk input distance = log2 rowsupp(K q_g)."""
import sys, random, numpy as np
np.seterr(all='ignore')
sys.path.insert(0, __import__('os').path.dirname(__file__))
from groups import G, z7, z2cube
from collections import Counter
grp, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
els = z7()[0] if grp == 'Z7' else z2cube()[0]
W = len(els); rng = random.Random(seed)
cnt = Counter(); found = []
for it in range(N):
    zs = (rng.randrange(1, 16), rng.randrange(1, 16))
    L = [G.rand_engine(rng) for _ in range(3)]
    t = G.mul(L[0], G.D[zs[0]], L[1], G.D[zs[1]], L[2])
    if G.dist_in(t) != 2: continue
    ti = G.inv(t); S = 1; d = []
    for g in range(1, W):
        Kq = G.mul(G.K, ti, els[g], t)
        if G.is_mono_any(G.mul(Kq, G.Ki)): S += 1
        else: d.append(G.dist_in(Kq))
    cnt[(S, d.count(4))] += 1
    if S < W and 4 not in d:
        found.append((it, zs, sorted(d)))
        if len(found) == 1:
            np.save(__import__('os').path.join(__import__('os').path.dirname(__file__), f'cex_t_{grp}.npy'), t)
print(grp, 'distance-2 terms sampled:', sum(cnt.values()), ' (|S|, #junk at distance 4):', dict(sorted(cnt.items())))
print('R\' condition fails (S != H, every junk cell nearer than 4):', len(found), found[:5])
