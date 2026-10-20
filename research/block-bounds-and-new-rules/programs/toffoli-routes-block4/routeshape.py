"""routeshape.py m -- for every route through a representative of each AGL-orbit of first backward steps:
T_K = its first backward step (at the kernel end) is a translation matching;
T_B = its first move (at the start end) is a translation, read through the transpose symmetry K^T = K (the transposed
      route's first backward step from K is a translation matching);
comp = number of connected components of the row/column support graph of each intermediate element (depth 1..m-1);
       an affine route has 2^(m-j) components at depth j (direct sum of 2^(m-j) kernels of block size j).
Weighted by orbit sizes (all invariants are unchanged by the route symmetry Lambda = right multiplication)."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, time, collections
import numpy as np
from rt import Game, A_, B_

m = int(sys.argv[1])
g = Game(m); n = g.n
cen = json.load(open("census_m%d.json" % m)  )

def routes_from(h, lev):
    if lev == 0:
        yield [h]; return
    for kk, (M, c) in g.children(h).items():
        if g.continues(c, lev - 1):
            for tail in routes_from(c, lev - 1):
                yield [h] + tail

def ncomp(h):
    S = np.abs(h) > 1e-12
    parent = list(range(2 * n))
    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in range(n):
        for c in np.nonzero(S[r])[0]:
            a, b = f(r), f(n + int(c))
            if a != b:
                parent[a] = b
    return len(set(f(x) for x in range(2 * n)))

def first_move_matching_T(chain):
    """chain = [h_m=K, ..., h_0]. The route's first move from B: h_1 = P (a + b s) h_0 (P monomial).
    Transposed route: K = K^T = h_0^T (a+b s^T) P^T (...); its first backward step from K is the coset of
    (a + b s^T)^{-1} h_0^{-T} K: read the pairs of rows of K it mixes."""
    h0, h1 = chain[-1], chain[-2]
    G1 = g.snap(h1 @ np.conj(h0).T)        # = P (a + b s): each row has two nonzeros
    T = G1.T                                # (a + b s^T) P^T
    # rows of K mixed by (T)^{-1}-type step: the transposed route's step (a + b s^T) acts on rows of h_0^T ... use
    # the matrix Y = (h1^T)^{-1} K = conj(h1) K  [h1 unitary], whose coset is the transposed route's first back step
    Y = g.snap(np.conj(h1) @ g.K)
    Z = g.snap(Y @ np.conj(g.K).T)          # Y K^-1: rows show which rows of K are combined
    pairs = set()
    for r in range(n):
        nz = tuple(np.nonzero(np.abs(Z[r]) > 1e-9)[0])
        assert len(nz) == 2, nz
        pairs.add(nz)
    assert len(pairs) == n // 2
    diffs = set(p ^ q for (p, q) in pairs)
    # sanity: Y really is one halving step back from K (rows of support n/2)
    assert (np.sum(np.abs(Y) > 1e-9, axis=1) == n // 2).all()
    return len(diffs) == 1

t0 = time.time()
dist = collections.Counter()
for r in cen:
    if not r["cont"]:
        continue
    M = [tuple(x) for x in r["rep"]]
    E0 = g.compat(g.K)
    h3 = g.back_step(g.K, M, E0)
    cnt = 0
    for tail in routes_from(h3, m - 1):
        chain = [g.K] + tail
        TB = first_move_matching_T(chain)
        comps = tuple(ncomp(chain[m - j]) for j in range(1, m))   # depth j element is chain[m-j]
        dist[(r["translation"], TB, comps)] += r["size"]
        cnt += 1
    assert cnt == r["nroutes"]
    print("rep size %d done (%d routes) t=%.1f" % (r["size"], cnt, time.time() - t0), flush=True)
print("(T_K, T_B, components at depth 1..m-1) : number of routes (all, weighted)")
for k in sorted(dist, key=str):
    print(k, dist[k])
print("total", sum(dist.values()))
tk = sum(v for k, v in dist.items() if k[0]); tb = sum(v for k, v in dist.items() if k[1])
both = sum(v for k, v in dist.items() if k[0] and k[1]); neither = sum(v for k, v in dist.items() if not k[0] and not k[1])
print("T_K:", tk, "T_B:", tb, "both:", both, "neither:", neither)
