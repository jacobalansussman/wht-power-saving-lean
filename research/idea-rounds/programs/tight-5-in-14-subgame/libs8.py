# libs8.py -- the singular subgame: G' = O(q0) = O6+(2) ~ S8 acting on 8 "points"; B' = Stab(L_Z) ~ AGL(3,2);
# h = t_(e1|e1) (a reflection = a transposition), K = h1 h2 h3; Young's orthogonal form for every irreducible of S8.
import itertools, numpy as np
from functools import lru_cache
from lib514 import *

NS = [v for v in range(1, NV) if q0(v) == 1]           # 28 nonsingular vectors <-> transpositions
assert len(NS) == 28
# points: maximal cliques of 7 pairwise non-orthogonal nonsingular vectors
cliques = set()
for v in NS:
    nb = [u for u in NS if u != v and BTAB[u, v] == 1]
    # the two cliques through v: v + (vectors u with B(u,v)=1) split into two classes of 6 by mutual non-orthogonality
    for u in nb:
        cl = frozenset([v, u] + [w for w in nb if w != u and w != (u ^ v) and BTAB[w, u] == 1])
        if len(cl) == 7 and all(BTAB[a, b] == 1 for a in cl for b in cl if a != b): cliques.add(cl)
POINTS = sorted(cliques, key=lambda c: sorted(c)); assert len(POINTS) == 8, len(POINTS)
PIDX = {c: i for i, c in enumerate(POINTS)}
def to_perm(g):
    """element of O(q0) (64-perm) -> permutation of the 8 points (array p with p[i] = image of point i)"""
    return np.array([PIDX[frozenset(int(g[v]) for v in c)] for c in POINTS], dtype=np.int64)
def pair_of(v):
    """nonsingular vector -> the transposition {i,j} (points whose cliques contain v)"""
    return tuple(i for i, c in enumerate(POINTS) if v in c)

def perm_word(p):
    """word in adjacent transpositions s_k = (k, k+1) (0-based k) with product = p (p[i] = image of i);
    returns list of k, to be multiplied left to right as matrices rho(s_k1) rho(s_k2) ... = rho(p)"""
    # bubble sort on the one-line notation of p^-1 ... use: p = s_{k1} s_{k2} ... s_{kr}
    q = list(range(8)); target = list(p); word = []
    # we build p by right-multiplying: cur = s_{k1}...s_{kj}; cur[i] = image of i
    # simpler: decompose via insertion: find word w with w applied to identity arrangement gives p
    arr = list(p)   # arr[i] = p(i); sort arr by adjacent swaps of POSITIONS; each swap = right mult by s_k
    w = []
    a = arr[:]
    for i in range(8):
        for j in range(7 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; w.append(j)
    # a = p o s_{w1} o s_{w2} ... = id  => p = s_{wr} ... s_{w1}  (each s is an involution)
    return w[::-1]

def compose(p, q): return p[q]          # (p o q)(i) = p(q(i))

# ---------------- Young's orthogonal form
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest
def syt(shape):
    """standard Young tableaux of shape, as dict entry -> (row, col)"""
    n = sum(shape); res = []
    def rec(k, filled, pos):
        if k > n: res.append(dict(pos)); return
        for r in range(len(shape)):
            c = filled[r]
            if c < shape[r] and (r == 0 or filled[r - 1] > c):
                filled[r] += 1; pos[k] = (r, c); rec(k + 1, filled, pos); del pos[k]; filled[r] -= 1
    rec(1, [0] * len(shape), {}); return res
class YoungRep:
    def __init__(self, shape):
        self.shape = shape; T = syt(shape); self.dim = len(T)
        key = lambda t: tuple(t[k] for k in range(1, 9))
        self.T = T; idx = {key(t): i for i, t in enumerate(T)}
        self.S = []
        for k in range(1, 8):            # s_k swaps k, k+1 (1-based entries) ; point indices 0..7 <-> entries 1..8
            M = np.zeros((self.dim, self.dim))
            for i, t in enumerate(T):
                (r1, c1), (r2, c2) = t[k], t[k + 1]
                if r1 == r2: M[i, i] = 1.0
                elif c1 == c2: M[i, i] = -1.0
                else:
                    a = (c2 - r2) - (c1 - r1)          # axial distance content(k+1) - content(k)
                    M[i, i] = 1.0 / a
                    t2 = dict(t); t2[k], t2[k + 1] = t[k + 1], t[k]; j = idx[key(t2)]
                    M[i, j] = np.sqrt(1 - 1.0 / a ** 2)
            self.S.append(M)
        self.cache = {}
    def __call__(self, p):
        kp = tuple(int(x) for x in p)
        if kp not in self.cache:
            M = np.eye(self.dim)
            for k in perm_word(p): M = M @ self.S[k]
            self.cache[kp] = M
        return self.cache[kp]
