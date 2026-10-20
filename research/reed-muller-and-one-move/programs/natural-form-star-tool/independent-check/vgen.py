"""vgen.py (v-rm2-exact): block-m version (m = 3, 4, 5) of vf.py + vrec.py, own code: field by carry-less
multiplication, own self-dual basis search, canonical targets mult(alpha^2k), forms, flats, star recursion with
pooling, clique, |Y|-W test, matching.  usage: vgen.py m < lines "k1,k2,..." [--allzero]"""
import sys, itertools, numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
POLYS = {3: 0b1011, 4: 0b10011, 5: 0b100101}
class Block:
    def __init__(self, m):
        self.m = m; q = 1 << m; n = q - 1; self.n = n; P = POLYS[m]
        def gmul(a, b):
            r = 0
            while b:
                if b & 1: r ^= a
                b >>= 1; a <<= 1
                if a & q: a ^= P
            return r
        def tr(x):
            s, y = 0, x
            for _ in range(m): s ^= y; y = gmul(y, y)
            return s
        self.gmul, self.tr = gmul, tr
        al = [1]
        for _ in range(n - 1): al.append(gmul(al[-1], 2))
        assert len(set(al)) == n
        # self-dual basis by exhaustive search (smallest lexicographic)
        cand = [x for x in range(1, q) if tr(gmul(x, x)) == 1]
        B = None
        for c in itertools.combinations(cand, m):
            if all(tr(gmul(c[i], c[j])) == 0 for i in range(m) for j in range(i + 1, m)):
                B = c; break
        assert B; self.B = B
        self.PAIRS = [(j, k) for j in range(m) for k in range(j, m)]; NF = 1 << len(self.PAIRS); self.NF = NF
        Z = np.arange(1, q)
        VALS = np.zeros((NF, n), dtype=np.int64)
        for p, (j, k) in enumerate(self.PAIRS):
            VALS ^= (((np.arange(NF) >> p) & 1)[:, None]) * (((Z >> j) & (Z >> k) & 1)[None, :])
        self.MASK = (VALS << np.arange(n)[None, :]).sum(axis=1); self.WT = VALS.sum(axis=1)
        QB = ((np.arange(NF)[:, None] >> np.arange(len(self.PAIRS))[None, :]) & 1)
        self.PAR = np.zeros((NF, n), dtype=np.int64)
        for k in range(n):
            beta = al[(2 * k) % n]
            T = [[tr(gmul(beta, gmul(B[i], B[j]))) for j in range(m)] for i in range(m)]
            tv = np.array([T[j][kk] for (j, kk) in self.PAIRS])
            self.PAR[:, k] = (QB @ tv) & 1
        PAT = (self.PAR << np.arange(n)[None, :]).sum(axis=1)
        pats = sorted(set(PAT.tolist())); cl = {p: i for i, p in enumerate(pats)}
        self.VID = np.array([cl[p] for p in PAT.tolist()]); self.NV = len(pats)
        self.CLSPAR = np.array([[(p >> k) & 1 for k in range(n)] for p in pats])
        def span(g):
            S = {0}
            for x in g: S |= {y ^ x for y in S}
            return S
        self.FL = []
        for d in range(m):
            subs = set()
            for g in itertools.combinations(range(1, q), d):
                S = span(g)
                if len(S) == 1 << d: subs.add(sum(1 << x for x in S))
            fl = set()
            for V in subs:
                pts = [x for x in range(q) if V >> x & 1]
                for a in range(q): fl.add(sum(1 << (x ^ a) for x in pts))
            self.FL.append(sorted(fl))
        FI = [{A: i for i, A in enumerate(f)} for f in self.FL]
        self.UP = []
        for d in range(m - 1):
            g = (1 << (m - d)) - 1; rows = []
            for A in self.FL[d]:
                pts = [x for x in range(q) if A >> x & 1]; Bs = set()
                for p in range(q):
                    if not A >> p & 1: Bs.add(A | sum(1 << (x ^ p ^ pts[0]) for x in pts))
                assert len(Bs) == g; rows.append([FI[d + 1][b] for b in Bs])
            self.UP.append(np.array(rows))
        self.DM = [np.array([(A >> 1) & ((1 << n) - 1) for A in f]) for f in self.FL]
        self._pre = {}
    def pre(self, F):
        if F in self._pre: return self._pre[F]
        n = self.n; A = int(self.MASK[F]); dirs = [z for z in range(1, n + 1) if A >> (z - 1) & 1]; d = len(dirs)
        G = np.arange(1, self.NF); G = G[G != F]; u = np.zeros(len(G), dtype=np.int64)
        for t, z in enumerate(dirs): u |= ((self.MASK[G] >> (z - 1)) & 1) << t
        T = np.full((1 << d, self.NV), 255, dtype=np.int64); np.minimum.at(T, (u, self.VID[G]), self.WT[G])
        for b in range(d):
            Tv = T.reshape(-1, 2, 1 << b, self.NV); np.minimum(Tv[:, 0], Tv[:, 1], out=Tv[:, 0])
        out = []
        for dd in range(self.m):
            S = self.DM[dd] & A; idx = np.zeros(len(S), dtype=np.int64)
            for t, z in enumerate(dirs): idx |= ((S >> (z - 1)) & 1) << t
            sz = np.array([bin(int(x)).count('1') for x in S])
            out.append((T[idx], sz, S))
        self._pre[F] = out; return out
def matching(D, Z):
    D = list(D); ix = {z: i for i, z in enumerate(D)}; r, c = [], []
    for x in D:
        for y in D:
            if (x, y) not in Z: r.append(ix[x]); c.append(ix[y])
    if not r: return 0
    M = csr_matrix((np.ones(len(r)), (r, c)), shape=(len(D), len(D)))
    return int((maximum_bipartite_matching(M, perm_type='column') >= 0).sum())
def maxclique(V, adj):
    best = []
    def bk(R, P):
        nonlocal best
        if len(R) > len(best): best = R[:]
        for v in sorted(P):
            if len(R) + len(P) <= len(best): return
            bk(R + [v], P & adj[v]); P = P - {v}
    bk([], set(V)); return best
def run(bl, ks, maxd=99, rounds=3, allzero=True):
    n, m = bl.n, bl.m; W = len(ks)
    w = bl.CLSPAR[:, ks].sum(axis=1); r = w[bl.VID]; s = bl.WT - r
    if (s[1:] < 0).any(): return ("D",)
    tests = [F for F in range(1, bl.NF) if r[F] == 0 and bl.WT[F] <= maxd and (allzero or bl.VID[F] == bl.VID[0])]
    Z = set()
    for rnd in range(1, rounds + 1):
        newz = False
        for F in tests:
            P = bl.pre(F); u = [np.minimum(P[dd][1], (P[dd][0] - w[None, :]).min(axis=1)) for dd in range(m)]
            if rnd > 1:
                for dd in range(m):
                    for i in np.nonzero(u[dd] > 0)[0]:
                        Sd = [z for z in range(1, n + 1) if int(P[dd][2][i]) >> (z - 1) & 1]
                        if any((x, y) in Z for x in Sd for y in Sd): u[dd][i] = min(u[dd][i], matching(Sd, Z))
            for dd in range(m - 2, -1, -1):
                g = (1 << (m - dd)) - 1; sz = P[dd][1]; Bu = u[dd + 1][bl.UP[dd]]; val = Bu.sum(axis=1) - 2 * sz
                if ((val < 0) & (sz > 0)).any(): return ("Z", rnd, F, dd)
                u[dd] = np.where(sz > 0, np.minimum(u[dd], np.minimum(val // (g - 2), Bu.min(axis=1))), u[dd])
            for dd in range(m):
                for i in np.nonzero((u[dd] == 0) & (P[dd][1] > 0))[0]:
                    Sd = [z for z in range(1, n + 1) if int(P[dd][2][i]) >> (z - 1) & 1]
                    for x in Sd:
                        for y in Sd:
                            if (x, y) not in Z: Z.add((x, y)); newz = True
        V = [x for x in range(1, n + 1) if (x, x) in Z]
        adj = {x: {y for y in V if y != x and (x, y) in Z and (y, x) in Z} for x in V}
        Y = maxclique(V, adj)
        if len(Y) > W: return ("C", rnd, len(Y))
        if matching(range(1, n + 1), Z) < n - W: return ("A", rnd)
        for Q in range(1, bl.NF):
            Dq = [z for z in range(1, n + 1) if int(bl.MASK[Q]) >> (z - 1) & 1]
            if any((x, y) in Z for x in Dq for y in Dq) and matching(Dq, Z) < s[Q]: return ("M", rnd, Q)
        if not newz: break
    global LASTZ; LASTZ = Z
    return ("P", len(Z))
if __name__ == "__main__":
    m = int(sys.argv[1]); bl = Block(m); print("m", m, "basis", bl.B, "flats", [len(f) for f in bl.FL], flush=True)
    for line in sys.stdin:
        t = line.split()
        if not t: continue
        ks = [int(x) for x in t[0].split(',')]
        print(t[0], " ".join(t[1:]), run(bl, ks), flush=True)
