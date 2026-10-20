"""sp8.py -- w5-blueprint: the frame-sum relaxation of block-4 programs lives on the dual polar graph of Sp(8,2).

Cells B g (B = engine letters, which contain the Paulis and the phases) <-> Lagrangians Lam_g = sigma(g)^-1 L_Z of
F_2^8, sigma = the symplectic image (conjugation of Paulis).  Input distance = 4 - dim(Lam cap L_Z).  One paid move with a
free letter takes a cell to ANY adjacent Lagrangian (dim of intersection 3) and never to itself.
Vectors: ints 0..255 (bits 0-3 x-part, 4-7 z-part).  Matrices: tuple of 8 column images.  Lagrangians: 256-bit masks.
16x16 matrices (D, K, engine letters, the Z7 / Z2^3 lifts) come from v-w3-blueprint/code/grp.py, groups.py (read-only use).
"""
import sys, os, itertools
import numpy as np
np.seterr(all='ignore')
_VC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vlib')   # unchanged copies of v-w3-blueprint/code/{grp,groups}.py
import importlib.util as _ilu
def _load(name):
    spec = _ilu.spec_from_file_location(name, os.path.join(_VC, name + '.py'))
    mod = _ilu.module_from_spec(spec); sys.modules[name] = mod; spec.loader.exec_module(mod); return mod
_load('grp'); _groups = _load('groups')
G, z7, z2cube = _groups.G, _groups.z7, _groups.z2cube

n, m = 16, 4
POP = [bin(i).count('1') for i in range(256)]
_POPT = np.array(POP, dtype=np.int64)

def omega(u, v):
    return POP[((u & 15) & (v >> 4)) ^ ((u >> 4) & (v & 15))] & 1

# ---------------------------------------------------------------- symplectic image of a 16x16 Clifford matrix
def pauli(a, b):
    M = np.zeros((n, n), dtype=complex)
    for y in range(n):
        M[y ^ a, y] = (-1) ** POP[b & y]
    return M

def ident_pauli(P):
    a = int(np.argmax(np.abs(P[:, 0])))
    base = P[a, 0]
    b = 0
    for j in range(m):
        y = 1 << j
        v = P[y ^ a, y] / base
        if abs(v + 1) < 1e-6:
            b |= 1 << j
        elif abs(v - 1) > 1e-6:
            return None
    Q = pauli(a, b) * base
    return (a, b) if np.allclose(P, Q, atol=1e-6) else None

def sigma(g):
    gi = np.conj(g).T
    cols = []
    for j in range(8):
        a, b = (1 << j, 0) if j < 4 else (0, 1 << (j - 4))
        r = ident_pauli(g @ pauli(a, b) @ gi)
        assert r is not None, 'not Clifford'
        cols.append(r[0] | (r[1] << 4))
    return tuple(cols)

# ---------------------------------------------------------------- F_2 linear algebra on 8-bit vectors
def app(S, v):
    out = 0
    j = 0
    while v:
        if v & 1:
            out ^= S[j]
        v >>= 1; j += 1
    return out

def table(S):
    lo = [0] * 16; hi = [0] * 16
    for v in range(16):
        lo[v] = app(S, v); hi[v] = app(S, v << 4)
    return [lo[v & 15] ^ hi[v >> 4] for v in range(256)]

def mul(S, T):          # (S T)(v) = S(T(v))
    return tuple(app(S, T[j]) for j in range(8))

def inv(S):
    img = {}
    for v in range(256):
        img[app(S, v)] = v
    return tuple(img[1 << j] for j in range(8))

IDM = tuple(1 << j for j in range(8))

def is_symp(S):
    for i in range(8):
        for j in range(8):
            if omega(S[i], S[j]) != omega(1 << i, 1 << j):
                return False
    return True

def span_mask(vs):
    s = [0]
    for v in vs:
        if not any(x == v for x in s):
            s = s + [x ^ v for x in s]
    mk = 0
    for x in s:
        mk |= 1 << x
    return mk

def mask_elems(mk):
    out = []
    v = 0
    while mk:
        if mk & 1:
            out.append(v)
        mk >>= 1; v += 1
    return out

def img_mask(T, mk):     # T = 256-table
    out = 0
    v = 0
    while mk:
        if mk & 1:
            out |= 1 << T[v]
        mk >>= 1; v += 1
    return out

def dimcap(a, b):
    c = bin(a & b).count('1')
    return c.bit_length() - 1

LZ = span_mask([1 << (4 + j) for j in range(4)])

# ---------------------------------------------------------------- all 2295 Lagrangians, adjacency
class Polar:
    def __init__(self):
        rng = __import__('random').Random(1)
        gens = [table(sigma(G.D[z])) for z in range(1, 16)] + [table(sigma(G.rand_engine(rng))) for _ in range(30)]
        seen = {LZ: 0}; L = [LZ]; fr = [LZ]
        while fr:
            nf = []
            for x in fr:
                for T in gens:
                    y = img_mask(T, x)
                    if y not in seen:
                        seen[y] = len(L); L.append(y); nf.append(y)
            fr = nf
        self.L = L; self.idx = seen
        assert len(L) == 2295, len(L)
        self.dZ = [4 - dimcap(x, LZ) for x in L]
        # adjacency: Lagrangians sharing a 3-space
        self.nb = [[] for _ in L]
        bycap = {}
        for i, x in enumerate(L):
            el = mask_elems(x)
            # 3-subspaces of x: kernels of the 15 nonzero functionals on a basis
            basis = basis_of(el)
            for f in range(1, 16):
                sub = [v for v in el if POP[f & coords(basis, v)] % 2 == 0]
                mk = 0
                for v in sub:
                    mk |= 1 << v
                bycap.setdefault(mk, []).append(i)
        for mk, lst in bycap.items():
            assert len(lst) == 3, len(lst)
            for a in lst:
                for b in lst:
                    if a != b:
                        self.nb[a].append(b)
        assert all(len(v) == 30 for v in self.nb)
        self.B = [basis_of(mask_elems(x)) for x in L]
        words = np.zeros((len(L), 4), dtype=np.uint64)
        for i, x in enumerate(L):
            for w in range(4):
                words[i, w] = (x >> (64 * w)) & ((1 << 64) - 1)
        cnt = np.zeros((len(L), len(L)), dtype=np.int64)
        for w in range(4):
            a = words[:, w][:, None] & words[:, w][None, :]
            # popcount of uint64 via bytes
            b = a.view(np.uint8).reshape(len(L), len(L), 8)
            cnt += _POPT[b].sum(axis=2, dtype=np.int64)
        self.D = (4 - np.log2(cnt).round().astype(int)).astype(np.int8)
    def perm(self, S):             # S: symplectic matrix (8 columns); returns the permutation of Lagrangian indices
        T = table(S)
        out = np.empty(len(self.L), dtype=np.int64)
        for i, b in enumerate(self.B):
            out[i] = self.idx[span_mask([T[v] for v in b])]
        return out
    def act(self, T, i):           # T: 256-table of a symplectic matrix; returns index of T(Lam_i)
        return self.idx[img_mask(T, self.L[i])]
    def d(self, i, j):
        return int(self.D[i, j])

def basis_of(el):
    b = []
    sp = {0}
    for v in el:
        if v not in sp:
            b.append(v); sp |= {x ^ v for x in sp}
    return b

def coords(basis, v):
    # coordinates of v in the basis (brute force, 4 vectors)
    for c in range(16):
        w = 0
        for j in range(len(basis)):
            if (c >> j) & 1:
                w ^= basis[j]
        if w == v:
            return c
    raise ValueError

# ---------------------------------------------------------------- the symmetry groups (symplectic images)
def groups():
    Kc = G.K; Ki = np.conj(G.K).T
    out = {}
    for name, els in (('Z7', z7()[0]), ('Z2^3', z2cube()[0])):
        c = [sigma(e) for e in els]
        cp = [sigma(G.mul(Ki, np.conj(e).T, Kc)) for e in els]    # c'_h = K^-1 c_h^-1 K
        out[name] = (c, cp)
    return out

def hindex_table(name, c):
    """group law: for Z7 els[h] = c^h; for Z2^3 els indexed by 3 bits (XOR)."""
    W = len(c)
    if name == 'Z7':
        add = lambda a, b: (a + b) % 7
        neg = lambda a: (-a) % 7
    else:
        add = lambda a, b: a ^ b
        neg = lambda a: a
    return W, add, neg

def characters(name):
    if name == 'Z7':
        w = np.exp(2j * np.pi / 7)
        return [[w ** (k * h) for h in range(7)] for k in range(7)]
    chars = []
    for k in range(8):
        chars.append([(-1) ** POP[k & h] for h in range(8)])
    return chars
