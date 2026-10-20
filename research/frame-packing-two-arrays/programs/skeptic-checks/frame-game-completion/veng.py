"""veng.py (v-pk2-complete, skeptic): my own engine at block m, written from the RULES (theory-guide RULES.md s.0)
and pk2-complete PROOF.md statements only; no import of pk2/fp-build code.

Group elements: n x n complex matrices (entries in Z[i,1/2]; float arithmetic is exact for them), keyed exactly.
Exact coefficients: Gaussian rationals as pairs of Fractions (class GQ).
Labels u = a | (b << m), P_u = X^a Z^b, X^a|x> = |x^a>, Z^b|x> = (-1)^(b.x)|x>.
Symplectic image gbar: 2m x 2m F2 matrix, column j = label of g P_(e_j) g^-1 (bit i < m: a_i, bit m+i: b_i).
"""
import numpy as np, itertools
from fractions import Fraction as Fr

def pc(x): return bin(x).count("1")

class GQ:
    __slots__ = ("r", "i")
    def __init__(s, r, i=Fr(0)): s.r = Fr(r); s.i = Fr(i)
    @staticmethod
    def of(z):
        z = complex(z); return GQ(Fr(z.real), Fr(z.imag))        # Fraction(float) is exact (dyadic)
    def __add__(s, o): return GQ(s.r + o.r, s.i + o.i)
    def __mul__(s, o): return GQ(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    def scale(s, q): return GQ(s.r * q, s.i * q)
    def iszero(s): return s.r == 0 and s.i == 0
    def c(s): return complex(float(s.r), float(s.i))
    def __repr__(s): return f"({s.r}+{s.i}i)"

class VEng:
    def __init__(s, m):
        s.m = m; n = s.n = 1 << m; s.NL = 4 ** m; s.mask = n - 1
        s.I = np.eye(n, dtype=complex)
        s.X = [s._X(a) for a in range(n)]; s.Z = [s._Z(b) for b in range(n)]
        s.P = [s.X[u & s.mask] @ s.Z[u >> m] for u in range(s.NL)]
        s.om = np.array([[(pc((w & s.mask) & (u >> m)) + pc((w >> m) & (u & s.mask))) & 1 for u in range(s.NL)]
                         for w in range(s.NL)], dtype=np.int64)
        s.sgn = 1 - 2 * s.om                                   # (-1)^omega(w,u)
        s.K = s.I.copy()
        for j in range(m): s.K = s.dirm(1 << j) @ s.K
        s.store = {}
    def _X(s, a):
        M = np.zeros((s.n, s.n), complex)
        for x in range(s.n): M[x ^ a, x] = 1
        return M
    def _Z(s, b): return np.diag([(-1.0) ** (pc(b & x) & 1) for x in range(s.n)]).astype(complex)
    def dirm(s, z): return (1 + 1j) / 2 * s.I + (1 - 1j) / 2 * s._X(z)
    # ---- exact keys
    def key(s, M):
        A = M * (1 << (s.m + 2)); R = np.round(A.real); Im = np.round(A.imag)
        assert np.abs(A.real - R).max() < 1e-9 and np.abs(A.imag - Im).max() < 1e-9, "not dyadic"
        k = (R.astype(np.int64).tobytes() + Im.astype(np.int64).tobytes())
        if k not in s.store: s.store[k] = M
        return k
    # ---- Pauli identification (own method: brute force against all labels)
    def ident(s, M):
        """M = c X^a Z^b: X^a Z^b|0> = |a>, then X^a M = c Z^b is diagonal; verified against the full matrix"""
        a = int(np.argmax(np.abs(M[:, 0])))
        Dm = s.X[a] @ M; c = Dm[0, 0]
        b = 0
        for i in range(s.m):
            if abs(Dm[1 << i, 1 << i] / c + 1) < 1e-9: b |= 1 << i
        u = a | (b << s.m)
        assert abs(abs(c) - 1) < 1e-9 and np.abs(M - c * s.P[u]).max() < 1e-9, "not a Pauli multiple"
        return c, u
    def gbar(s, g):
        gi = g.conj().T
        Sb = np.zeros((2 * s.m, 2 * s.m), dtype=np.int64)
        for j in range(2 * s.m):
            _, u = s.ident(g @ s.P[1 << j] @ gi)
            for i in range(2 * s.m): Sb[i, j] = (u >> i) & 1
        return Sb
    def act(s, Sb, u):
        v = np.array([(u >> i) & 1 for i in range(2 * s.m)]); w = Sb @ v % 2
        return int(sum(int(w[i]) << i for i in range(2 * s.m)))
    @staticmethod
    def skey(Sb): return Sb.astype(np.uint8).tobytes()
    def in_P(s, Sb): return not Sb[:s.m, s.m:].any()       # a-part of images of (0|b) vanishes
    # ---- cocycles and D(g)
    def mu4(s, c):
        for j in range(4):
            if abs(c - 1j ** j) < 1e-9: return j
        raise AssertionError("not in mu4: %r" % c)
    def Dk(s, g):
        """coefficient vectors over [P_u] of d_k(g), k = 0..3 (PROOF Theorem 1 construction)"""
        gi = g.conj().T                                       # g unitary: exact inverse
        d0 = np.zeros(s.NL, complex); d2 = np.zeros(s.NL, complex)
        for w in range(s.NL):
            c0, _ = s.ident(g @ s.P[w] @ gi); s.mu4(c0); d0[w] = 1 / c0
            c2, _ = s.ident(g @ s.P[w] @ g.T); s.mu4(c2); d2[w] = 1 / c2
        D0 = (d0 @ s.sgn) / s.NL                              # sum_w d^-1 4^-m sum_u (-1)^om(w,u) [P_u]
        D2 = (d2 @ s.sgn) / s.NL
        D1 = np.array([np.trace(s.P[u].conj().T @ gi) / s.n for u in range(s.NL)])
        D3 = np.array([np.trace(s.P[u].conj().T @ g.T) / s.n for u in range(s.NL)])
        return [D0, D1, D2, D3]
    def D(s, g, mode='full'):
        """dict {(u, j): GQ}: the C[N] element sum_k d_k e_k; group element (u, j) = i^j P_u.
        mode 'f1': fp-build's f = 1 lift alone, Pexp(g^-1) with no block projection."""
        Ds = s.Dk(g); out = {}
        if mode == 'f1':
            for u in range(s.NL):
                if abs(Ds[1][u]) > 1e-12: out[(u, 0)] = GQ.of(Ds[1][u])
            return out
        for k in range(4):
            for u in range(s.NL):
                c = Ds[k][u]
                if abs(c) < 1e-12: continue
                cq = GQ.of(c)
                for j in range(4):
                    t = cq * GQ.of((1j) ** ((-j * k) % 4)); t = t.scale(Fr(1, 4))
                    out[(u, j)] = out.get((u, j), GQ(0)) + t
        return {k: v for k, v in out.items() if not v.iszero()}
    def nmat(s, u, j): return (1j ** j) * s.P[u]
    # ---- my lifts of P from legitimate free letters
    def phase_letter(s, z, c):
        return np.diag([1j ** ((pc(z & x) + c) & 1) for x in range(s.n)])
    def addr_mat(s, A):
        """A: list of m column images (ints); |x> -> |A x>"""
        M = np.zeros((s.n, s.n), complex)
        for x in range(s.n):
            y = 0
            for k in range(s.m):
                if (x >> k) & 1: y ^= A[k]
            M[y, x] = 1
        return M
    def build_Ptable(s):
        """every p in P -> a word in legitimate free letters (phase letters, address maps) and its matrix"""
        m = s.m; s.Ptab = {}
        gl = []
        for cols in itertools.product(range(1, s.n), repeat=m):
            # invertible?
            span = {0}
            for c in cols: span |= {x ^ c for x in span}
            if len(span) == s.n: gl.append(list(cols))
        pairs = [(i, j) for i in range(m) for j in range(i, m)]
        for A in gl:
            Ad = s.addr_mat(list(A))
            for bits in range(1 << len(pairs)):
                M = Ad.copy(); word = [('addr', tuple(A))]
                for t, (i, j) in enumerate(pairs):
                    if not (bits >> t) & 1: continue
                    if i == j:
                        L = [('ph', 1 << i, 0)]
                    else:   # CZ_ij = ph(e_i) ph(e_j) ph(e_i+e_j)^3
                        L = [('ph', 1 << i, 0), ('ph', 1 << j, 0)] + [('ph', (1 << i) | (1 << j), 0)] * 3
                    for (_, z, c) in L:
                        M = s.phase_letter(z, c) @ M; word.append(('ph', z, c))
                Sb = s.gbar(M); assert s.in_P(Sb)
                k = s.skey(Sb)
                if k not in s.Ptab: s.Ptab[k] = (M, word)
        return len(s.Ptab)
