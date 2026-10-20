# vgamma.py -- v-nc-gain's own frame graph Gamma_m (independent of nc-gain's dpg.py).
# Lagrangians of F_2^(2m) found by recursive extension of isotropic subspaces (canonical: reduced echelon basis);
# adjacency = intersection of dimension m-1; distances by BFS (then checked against m - dim(E cap F)).
# Vector encoding (a|b) -> a | (b << m) and form w = a.b' + b.a' are the SAME as dpg.py, only so that frame indices can
# be mapped onto nc-gain's by subspace equality (the construction itself shares no code).
import numpy as np, itertools, collections, sys

def w(u, v, m):
    M = (1 << m) - 1
    return (bin((u & M) & (v >> m)).count('1') + bin((u >> m) & (v & M)).count('1')) & 1

def rref_add(basis, v):
    """basis: dict pivot->vector (reduced). return new basis with v added, or None if v in span."""
    for p, b in basis.items():
        if v >> p & 1: v ^= b
    if v == 0: return None
    p = v.bit_length() - 1
    nb = {}
    for q, b in basis.items():
        nb[q] = b ^ v if b >> p & 1 else b
    nb[p] = v
    return nb

def lagrangians(m):
    N = 1 << (2 * m)
    seen = set(); out = []
    level = {frozenset(): {}}
    for k in range(m):
        nxt = {}
        for key, basis in level.items():
            vecs = list(basis.values())
            for v in range(1, N):
                if w(v, v, m): continue          # alternating form: always 0, kept for clarity
                if any(w(v, b, m) for b in vecs): continue
                nb = rref_add(basis, v)
                if nb is None: continue
                kk = frozenset(nb.values())
                if kk not in nxt: nxt[kk] = nb
        level = nxt
    for key, basis in level.items():
        vecs = list(basis.values()); S = {0}
        for v in vecs: S |= {x ^ v for x in S}
        out.append(frozenset(S))
    return out

class VG:
    def __init__(self, m):
        self.m = m
        L = lagrangians(m)
        # order: put the raw frame {(0|b)} first, rest sorted by min nonzero element pattern (arbitrary but fixed)
        L0 = frozenset(b << m for b in range(1 << m))
        L.sort(key=lambda S: (S != L0, sorted(S)))
        self.L = L; n = self.n = len(L)
        bits = 1 << (2 * m)
        Mx = np.zeros((n, bits), np.uint8)
        for i, S in enumerate(L): Mx[i, list(S)] = 1
        inter = Mx.astype(np.int32) @ Mx.T.astype(np.int32)        # |E cap F| = 2^dim
        dimI = np.log2(inter).round().astype(np.int16)
        self.dimI = dimI
        A = (dimI == m - 1)
        self.adj = [np.nonzero(A[i])[0] for i in range(n)]
        # BFS distances
        D = np.full((n, n), -1, np.int16)
        srcs = range(n) if n < 1000 else list(range(0, n, 97))
        for s in srcs:
            D[s, s] = 0; fr = [s]; d = 0
            while fr:
                d += 1; nf = []
                for x in fr:
                    for y in self.adj[x]:
                        if D[s, y] < 0: D[s, y] = d; nf.append(y)
                fr = nf
        if n >= 1000:
            Dc = (m - dimI).astype(np.int16)
            assert (D[srcs] == Dc[srcs]).all(), 'BFS sample disagrees with m - dim(cap)'
            D = Dc
        self.D = D
        self.edges = [(u, int(v)) for u in range(n) for v in self.adj[u] if u < v]

def check(m):
    g = VG(m)
    deg = collections.Counter(len(a) for a in g.adj)
    diam = int(g.D.max())
    ok_dist = bool((g.D == (m - g.dimI)).all())
    opp = int((g.D == m).sum())
    print(f"m={m}: frames {g.n}, degrees {dict(deg)}, edges {len(g.edges)}, diameter {diam}, ordered opposite pairs {opp}, "
          f"d = m - dim(cap) everywhere: {ok_dist}; distance profile from frame 0 {dict(collections.Counter(g.D[0].tolist()))}")
    return g

if __name__ == '__main__':
    for m in map(int, sys.argv[1:]): check(m)
