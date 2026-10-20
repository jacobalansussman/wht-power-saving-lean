"""junkdist.py -- where do the junk cells [K q_g] of a shared move's content term t sit?  q_g = t^-1 c_g t.
For t a raw-input letter b (shared START, Theorem R): K q_g is monomial-times-K... full row support always.
For t = D1 a (an O1 output: shared move AFTER the first free orbit): input distance of the junk cell = log2 rowsupp(K q_g)."""
import sys
from symsolve import *
from collections import Counter
W, H, n = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
R = Run(4, W, H, 1); G = R.G; rng = random.Random(2)
def rs(M):
    return int(round(np.log2(G.rowsupp(M).max())))
for kind in ('start b', 'O1 output D1 a'):
    cnt = Counter()
    for _ in range(n):
        a = rand_free_b(G, rng, 4)
        t = a if kind == 'start b' else G.mul(G.D1, a)
        ti = np.conj(t).T
        Kc = np.conj(G.K).T
        def lab(h):
            M = G.mul(G.K, ti, h[0], t)
            if G.is_monomial(G.mul(M, Kc)):
                return 'S'                       # g in S: the cell is the kernel cell itself
            return rs(M)
        dists = tuple(sorted((lab(h) for h in R.Hc[1:]), key=str))
        cnt[dists] += 1
    print(kind, ': input distances of the junk cells [K q_g], g != 0 (sorted per term; S = g in S, i.e. the kernel cell itself):', cnt.most_common(6), flush=True)
