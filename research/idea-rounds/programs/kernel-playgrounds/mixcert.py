# playgrounds/mixcert.py -- monomial-invariant coordinate certificates in V (x) S^k V and S^a V (x) Lambda^b V etc.
# A basis vector of S^a V (x) Lambda^b V is (multiset A, set B); Y = span of those whose "Venn data" (for each
# element: multiplicity in A, membership in B) satisfies a monomial-invariant condition (a predicate on the
# multiset of per-element patterns). Every such Y is kept by all monomial matrices. Theorem R: moves/array >= dem/cap.
# Usage: python3 -I mixcert.py
import sys, itertools
from collections import Counter
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge, kron_list
P = 2147483647

def sympow(g, k, p, msets):
    """S^k(g) mod p in the monomial basis (multisets of size k): column U -> coefficients of prod_{u in U} (g e_u)."""
    n = g.shape[0]; idx = {T: i for i, T in enumerate(msets)}
    out = np.zeros((len(msets), len(msets)), dtype=np.int64)
    for j, U in enumerate(msets):
        poly = {(): 1}
        for u in U:
            col = [(t, int(g[t, u]) % p) for t in range(n) if g[t, u] % p]
            new = {}
            for mono, c in poly.items():
                for t, gv in col:
                    key = tuple(sorted(mono + (t,)))
                    new[key] = (new.get(key, 0) + c * gv) % p
            poly = new
        for mono, c in poly.items():
            out[idx[mono], j] = (out[idx[mono], j] + c) % p
    return out

def certs(K, m, kind, a, b, p=P):
    n = 2 ** m
    K = np.array(K, dtype=object); I = np.eye(2, dtype=object)
    K1 = kron_list([K] + [I] * (m - 1)); T = kron_list([K] * m)
    MA = list(itertools.combinations_with_replacement(range(n), a)) if kind[0] == 'S' else list(itertools.combinations(range(n), a))
    MB = list(itertools.combinations_with_replacement(range(n), b)) if kind[1] == 'S' else list(itertools.combinations(range(n), b))
    def rep(g, which, k, M):
        return sympow(g, k, p, M) if which == 'S' else wedge(g, k, p, M)
    mats = {}
    for name, g in (('cap', K1), ('dem', T)):
        gp = np.array([[int(v) % p for v in row] for row in g], dtype=np.int64)
        mats[name] = np.kron(rep(gp, kind[0], a, MA), rep(gp, kind[1], b, MB)) % p
    # Venn signature: multiset over elements of (mult in A, mult in B), elements with (0,0) dropped
    sig = []
    for A in MA:
        ca = Counter(A)
        for B in MB:
            cb = Counter(B)
            s = Counter((ca[e], cb[e]) for e in set(ca) | set(cb))
            sig.append(tuple(sorted(s.items())))
    sigs = sorted(set(sig)); sig = np.array([sigs.index(s) for s in sig])
    return mats, sig, sigs

def scan(K, m, kind, a, b, maxsub=4096):
    mats, sig, sigs = certs(K, m, kind, a, b)
    best = (0, None)
    nsig = len(sigs)
    subsets = itertools.chain.from_iterable(itertools.combinations(range(nsig), r) for r in range(1, nsig))
    cnt = 0
    for S in subsets:
        cnt += 1
        if cnt > maxsub: break
        inY = np.isin(sig, S)
        cap = rank_mod(mats['cap'][np.ix_(~inY, inY)], P)
        if cap == 0: continue
        dem = rank_mod(mats['dem'][np.ix_(~inY, inY)], P)
        if dem / cap > best[0]: best = (dem / cap, [sigs[i] for i in S], cap, dem, int(inY.sum()))
    return best, nsig, len(sig)

if __name__ == '__main__':
    m = 3
    for kind, a, b in (('VS', 1, 2), ('VS', 1, 3), ('VS', 1, 4), ('SL', 2, 2), ('SL', 2, 3), ('SS', 2, 2)):
        for name, K in (('H', [[1, 1], [1, -1]]), ('L', [[1, 0], [1, 1]])):
            kk = ('L' if kind[0] == 'V' else kind[0], kind[1])
            best, nsig, dim = scan(K, m, kk, a, b)
            print(f"{name} {kind} a={a} b={b} dim={dim} #sig={nsig}: best ratio {best[0]:.4f} cap={best[2]} dem={best[3]} dimY={best[4]} Y={best[1]}", flush=True)
