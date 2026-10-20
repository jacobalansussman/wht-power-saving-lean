"""vgeo.py (v-w7-b5, skeptic's own code, pure F_2): dual polar geometry of Sp(10,2) and orbit statuses for the decoupled
2-in-8 of w7-b5.  Vectors are 10-bit ints: low 5 bits = e-coordinates u, high 5 bits = f-coordinates w;
omega((u,w),(u',w')) = u.w' + w.u'.  A frame (Lagrangian) is a frozenset of its 32 vectors.  d(X,Y) = 5 - log2|X & Y|.
M0 = span(e1..e5).  Stab(M0) = {(u,w) -> (A u + A S w, A^-T w) : A in GL5, S symmetric} (all of it).
usage: vgeo.py seed nconf nlist  (nlist random Stab(M0) elements per list; sigma found by birthday matching)"""
import sys, random

def omega(x, y):
    return bin(((x & 31) & (y >> 5)) ^ ((x >> 5) & (y & 31))).count('1') & 1

def span(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {a ^ v for a in s}
    return frozenset(s)

def basis(X):
    b = []; s = {0}
    for v in sorted(X):
        if v not in s:
            b.append(v); s |= {a ^ v for a in s}
    return b

def dist(X, Y):
    return 5 - (len(X & Y).bit_length() - 1)

def isotropic(X):
    b = basis(X)
    return all(omega(a, c) == 0 for a in b for c in b) and len(X) == 32

def perp(vs):
    return [v for v in range(1024) if all(omega(v, h) == 0 for h in vs)]

def neighbours(X):
    b = basis(X); out = set()
    for f in range(1, 32):
        # hyperplane H = {sum c_i b_i : sum c_i f_i = 0}
        Hv = []
        for c in range(32):
            if bin(c & f).count('1') % 2 == 0:
                v = 0
                for i in range(5):
                    if (c >> i) & 1:
                        v ^= b[i]
                Hv.append(v)
        Hb = basis(frozenset(Hv))
        assert len(Hb) == 4
        for v in perp(Hb):
            if v not in X:
                out.add(span(Hb + [v]))
    assert len(out) == 62, len(out)
    return out

# ---- group elements: list of images of the 10 basis vectors (bit i)
def apply(g, x):
    y = 0; i = 0
    while x:
        if x & 1:
            y ^= g[i]
        x >>= 1; i += 1
    return y

def gmul(g, h):  # g o h
    return [apply(g, h[i]) for i in range(10)]

def gframe(g, X):
    return span([apply(g, v) for v in basis(X)])

def mat5_inv(A):  # A as list of 5 column ints (5-bit); returns inverse or None
    # rows representation
    M = [[(A[j] >> i) & 1 for j in range(5)] + [1 if i == k else 0 for k in range(5)] for i in range(5)]
    r = 0
    for c in range(5):
        p = next((i for i in range(r, 5) if M[i][c]), None)
        if p is None:
            return None
        M[r], M[p] = M[p], M[r]
        for i in range(5):
            if i != r and M[i][c]:
                M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        r += 1
    inv_cols = []
    for j in range(5):
        col = 0
        for i in range(5):
            if M[i][5 + j]:
                col |= 1 << i
        inv_cols.append(col)
    return inv_cols

def m5apply(A, u):
    y = 0
    for i in range(5):
        if (u >> i) & 1:
            y ^= A[i]
    return y

def transpose5(A):
    return [sum((((A[j] >> i) & 1) << j) for j in range(5)) for i in range(5)]

def rand_stab_M0(rng):
    while True:
        A = [rng.randrange(32) for _ in range(5)]
        Ai = mat5_inv(A)
        if Ai is not None:
            break
    AiT = transpose5(Ai)
    S = [[0] * 5 for _ in range(5)]
    for i in range(5):
        for j in range(i, 5):
            S[i][j] = S[j][i] = rng.randrange(2)
    Scol = [sum(S[i][j] << i for i in range(5)) for j in range(5)]
    g = []
    for k in range(5):  # e_k -> (A e_k, 0)
        g.append(A[k])
    for k in range(5):  # f_k -> (A S f_k, A^-T f_k)
        g.append(m5apply(A, Scol[k]) | (m5apply(AiT, 1 << k) << 5))
    return g

def is_symplectic(g):
    return all(omega(g[i], g[j]) == omega(1 << i, 1 << j) for i in range(10) for j in range(10))

def ginv(g):
    # brute force: invert 10x10 over F2
    M = [[(g[j] >> i) & 1 for j in range(10)] + [1 if i == k else 0 for k in range(10)] for i in range(10)]
    r = 0
    for c in range(10):
        p = next(i for i in range(r, 10) if M[i][c])
        M[r], M[p] = M[p], M[r]
        for i in range(10):
            if i != r and M[i][c]:
                M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        r += 1
    return [sum(M[i][10 + j] << i for i in range(10)) for j in range(10)]

def order(g, cap=5000):
    I = [1 << i for i in range(10)]
    h = g; n = 1
    while h != I:
        h = gmul(g, h); n += 1
        if n > cap:
            return None
    return n
