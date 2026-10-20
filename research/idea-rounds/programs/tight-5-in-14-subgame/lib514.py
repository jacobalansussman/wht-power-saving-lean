# lib514.py -- block-3 frame game: Sp(6,2) on F_2^6, the parabolic P0 = Stab(L_Z), the Levi L = GL(3,2),
# the move h = t_(e1|e1), the target K = h1 h2 h3, and explicit irreducibles V_7, V_27, V_56 built from the
# 28 minus-type quadratic forms (own construction; only the character table is borrowed from c2-finite-group/gap).
# Vectors: v = x1 + 2 x2 + 4 x3 + 8 z1 + 16 z2 + 32 z3; B(u,v) = x_u.z_v + z_u.x_v.
# Group elements: int arrays g of length 64 with g[v] = image of v; composition comp(g,h) = g o h.
import itertools, json, warnings
import numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)  # spurious FP flags from macOS Accelerate matmul

NV = 64
def pc(x): return bin(x).count('1')
def Bf(u, v): return (pc((u & 7) & (v >> 3)) + pc((u >> 3) & (v & 7))) & 1
BTAB = np.array([[Bf(u, v) for v in range(NV)] for u in range(NV)], dtype=np.int64)
def trans(w): return np.array([v ^ (w if BTAB[v, w] else 0) for v in range(NV)], dtype=np.int64)
def comp(g, h): return g[h]
def inv(g):
    r = np.empty_like(g); r[g] = np.arange(NV); return r
def key(g): return g.astype(np.int8).tobytes()
ID = np.arange(NV, dtype=np.int64)
def gmat(g):
    """6x6 F_2 matrix (columns = images of basis vectors 1,2,4,8,16,32)"""
    return np.array([[(int(g[1 << j]) >> i) & 1 for j in range(6)] for i in range(6)], dtype=np.int64)
def rank2(M):
    M = (np.array(M, dtype=np.int64) % 2).copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r
def order(g):
    x = g.copy(); k = 1
    while not np.array_equal(x, ID): x = comp(g, x); k += 1
    return k
def span(vs):
    S = {0}
    for v in vs: S |= {x ^ v for x in S}
    return frozenset(S)

X = lambda i: 1 << (i - 1)
Z = lambda i: 1 << (2 + i)
EE = lambda i: X(i) | Z(i)
TRANS = {w: trans(w) for w in range(1, NV)}
h1 = TRANS[EE(1)]; h2 = TRANS[EE(2)]; h3 = TRANS[EE(3)]
H = h1                                   # the paid move (frame-game image of H_1)
K = comp(comp(h1, h2), h3)               # the target (frame-game image of H^{(x)3})
LZ = span([Z(1), Z(2), Z(3)]); LX = span([X(1), X(2), X(3)])
A_GENS = [TRANS[Z(1)], TRANS[Z(2)], TRANS[Z(3)], TRANS[Z(1) | Z(2)], TRANS[Z(1) | Z(3)], TRANS[Z(2) | Z(3)]]
def inv3(A):
    A = np.array(A) % 2; M = np.concatenate([A, np.eye(3, dtype=int)], axis=1)
    for c in range(3):
        p = next(r for r in range(c, 3) if M[r, c]); M[[c, p]] = M[[p, c]]
        for r in range(3):
            if r != c and M[r, c]: M[r] ^= M[c]
    return M[:, 3:]
def levi(Amat):
    """(x|z) -> (A x | A^{-T} z): the Levi element of the address map A (3x3 over F_2)"""
    Am = np.array(Amat) % 2; AinvT = inv3(Am).T
    def ap(Mt, x): return sum(((int(Mt[i] @ np.array([(x >> j) & 1 for j in range(3)])) & 1) << i) for i in range(3))
    return np.array([ap(Am, v & 7) | (ap(AinvT, v >> 3) << 3) for v in range(NV)], dtype=np.int64)
ALL3 = [np.array(m).reshape(3, 3) for m in itertools.product([0, 1], repeat=9)]
GL3 = [A for A in ALL3 if round(abs(np.linalg.det(A))) % 2 == 1 and rank2(A) == 3]
assert len(GL3) == 168
LEVI = [levi(A) for A in GL3]
def gen_closure(gens, limit=10**7):
    el = {key(ID): ID}; fr = [ID]
    while fr:
        nx = []
        for x in fr:
            for s in gens:
                y = comp(s, x); k = key(y)
                if k not in el: el[k] = y; nx.append(y)
        fr = nx
        if len(el) > limit: raise RuntimeError('too big')
    return list(el.values())

# ---------------- quadratic forms q_a(v) = x.z + B(a, v); minus type iff x_a . z_a = 1
def q0(v): return pc((v & 7) & (v >> 3)) & 1
QTAB = np.array([[(q0(v) + BTAB[a, v]) & 1 for v in range(NV)] for a in range(NV)], dtype=np.int64)  # QTAB[a, v]
QIDX = {QTAB[a].tobytes(): a for a in range(NV)}
MINUS = [a for a in range(NV) if q0(a) == 1]; PLUS = [a for a in range(NV) if q0(a) == 0]
assert len(MINUS) == 28 and len(PLUS) == 36
def form_act(g):
    """g . q_a = q_{b}: returns array b[a] (q o g^{-1})"""
    gi = inv(g)
    return np.array([QIDX[QTAB[a][gi].tobytes()] for a in range(NV)], dtype=np.int64)
MI = {a: i for i, a in enumerate(MINUS)}
a0 = MINUS[0]
# coset representatives gamma_a (gamma_a . q_a0 = q_a), BFS over transvections
GAM = {a0: ID}; order_ = [a0]
for a in order_:
    for w in range(1, NV):
        g = TRANS[w]; b = int(form_act(g)[a])
        if b not in GAM: GAM[b] = comp(g, GAM[a]); order_.append(b)
assert len(GAM) == 28
GAMI = {a: inv(g) for a, g in GAM.items()}
def dickson(g): return rank2((gmat(g) + np.eye(6, dtype=np.int64)) % 2) & 1
def signed28(g):
    """Ind_{O(q_a0)}^G (Dickson sign): monomial 28x28 matrix"""
    fa = form_act(g); S = np.zeros((28, 28))
    for a in MINUS:
        b = int(fa[a]); s = comp(GAMI[b], comp(g, GAM[a]))
        assert int(form_act(s)[a0]) == a0
        S[MI[b], MI[a]] = -1.0 if dickson(s) else 1.0
    return S
def perm28(g):
    fa = form_act(g); P = np.zeros((28, 28))
    for a in MINUS: P[MI[int(fa[a])], MI[a]] = 1.0
    return P

def eigbasis(Zs, lam, tol=1e-6):
    """orthonormal basis of ker(Zs - lam) for a real symmetric Zs"""
    w, V = np.linalg.eigh(Zs)
    return V[:, np.abs(w - lam) < tol], w

def build_irreps(verbose=False):
    Zs = sum(signed28(TRANS[w]) for w in range(1, NV)); B7, w7 = eigbasis(Zs, -45.0)
    Zp = sum(perm28(TRANS[w]) for w in range(1, NV)); B27, w27 = eigbasis(Zp, 63 * 15 / 27)
    assert B7.shape[1] == 7 and B27.shape[1] == 27, (B7.shape, B27.shape)
    def r7(g): return B7.T @ signed28(g) @ B7
    def r27(g): return B27.T @ perm28(g) @ B27
    Z189 = sum(np.kron(r7(TRANS[w]), r27(TRANS[w])) for w in range(1, NV))
    Z189 = (Z189 + Z189.T) / 2
    B56, w189 = eigbasis(Z189, -27.0)
    assert B56.shape[1] == 56, B56.shape
    if verbose:
        print('V_28signed transvection-sum eigenvalues', sorted(set(np.round(np.linalg.eigvalsh(Zs), 6))))
        print('V_7 (x) V_27 transvection-sum eigenvalues', sorted(set(np.round(w189, 6))))
    return B7, B27, B56

class Irreps:
    def __init__(self, verbose=False):
        self.B7, self.B27, self.B56 = build_irreps(verbose)
        self.cache = {}
    def r7(self, g): return self.B7.T @ signed28(g) @ self.B7
    def r27(self, g): return self.B27.T @ perm28(g) @ self.B27
    def r56(self, g):
        k = key(g)
        if k not in self.cache:
            self.cache[k] = self.B56.T @ np.kron(self.r7(g), self.r27(g)) @ self.B56
        return self.cache[k]

def rank(M, tol=1e-8):
    if M.size == 0: return 0
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol * max(1.0, s[0])).sum())
