# Skeptic v-ex-dimension-compute, 14 Oct 2026.  Independent mod-p implementation of the engine letters (typed from
# the definitions in proofs/ladder/ladder_check.py and pf-algebraic Lemma 6, NOT imported from dimcount/bsearch).
import numpy as np, itertools
def setup(p):
    assert p % 4 == 1
    for g in range(2, p):
        I = pow(g, (p - 1) // 4, p)
        if I * I % p == p - 1: break
    return I, pow(2, p - 2, p)
def perm(n, f):                 # (P v)(x) = v(f(x))
    M = np.zeros((n, n), np.int64)
    for x in range(n): M[x, f(x)] = 1
    return M
def letters(m, p):
    I, h = setup(p); n = 1 << m
    free = []
    for j in range(m): free.append(perm(n, lambda x, j=j: x ^ (1 << j)))                       # X-shifts
    for j in range(1, m):                                                                      # swaps of bit 0 and bit j
        def sw(x, j=j):
            b0, bj = x & 1, (x >> j) & 1
            return (x & ~(1 | (1 << j))) | (bj) | (b0 << j)
        free.append(perm(n, sw))
    if m >= 2: free.append(perm(n, lambda x: x ^ ((x & 1) << 1)))                           # transvection
    for j in range(m): free.append(np.diag([I if (x >> j) & 1 else 1 for x in range(n)]).astype(np.int64))  # S_j
    free.append((I * np.eye(n, dtype=np.int64)) % p)                                           # scalar i
    ca, cb = (1 + I) * h % p, (1 - I) * h % p
    def dirm(j): return (ca * np.eye(n, dtype=np.int64) + cb * perm(n, lambda x: x ^ (1 << j))) % p
    Wr = [np.eye(n, dtype=np.int64)]
    for j in range(m): Wr.append(dirm(j) @ Wr[-1] % p)
    return free, Wr
def mm(a, b, p): return np.matmul(a, b) % p
def ext2(g, p):
    n = g.shape[0]; S = list(itertools.combinations(range(n), 2)); k = len(S)
    M = np.zeros((k, k), np.int64)
    for a, (i, j) in enumerate(S):
        for b, (c, d) in enumerate(S):
            M[a, b] = (g[i, c] * g[j, d] - g[i, d] * g[j, c]) % p
    return M
def rep(name, g, p):
    if name == 'V': return g % p
    if name == 'VV': return np.kron(g, g) % p
    if name == 'VL2': return np.kron(g, ext2(g, p)) % p
    raise ValueError
def rank_mod(M, p):
    M = M.copy() % p; r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if nz.size == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = M[r] * pow(int(M[r, c]), p - 2, p) % p
        col = M[r + 1:, c].copy(); nzr = np.nonzero(col)[0]
        if nzr.size:
            M[r + 1 + nzr] = (M[r + 1 + nzr] - np.outer(col[nzr], M[r]) % p) % p
        r += 1
    return r
class Echelon:
    """incremental row-echelon basis mod p (rows reduced, pivots normalised to 1)"""
    def __init__(self, ncols, p):
        self.p = p; self.rows = np.zeros((0, ncols), np.int64); self.piv = []
    def reduce(self, X):
        p = self.p; X = X % p
        for k, c in enumerate(self.piv):
            f = X[:, c].copy(); nz = np.nonzero(f)[0]
            if nz.size: X[nz] = (X[nz] - np.outer(f[nz], self.rows[k]) % p) % p
        return X
    def add(self, X):
        """add rows of X; returns the new independent rows (reduced), in echelon form"""
        p = self.p; X = self.reduce(X); new = []
        for i in range(X.shape[0]):
            v = X[i]
            for k in range(len(new)):
                c = new[k][1]
                if v[c]: v = (v - v[c] * new[k][0]) % p
            nz = np.nonzero(v)[0]
            if nz.size == 0: continue
            c = nz[0]; v = v * pow(int(v[c]), p - 2, p) % p
            # also clear this pivot from the later rows of X (lazy: done in the loop above)
            new.append((v, c))
        if new:
            self.rows = np.vstack([self.rows] + [v[None] for v, c in new]); self.piv += [c for v, c in new]
        return [v for v, c in new]
    def dim(self): return len(self.piv)
