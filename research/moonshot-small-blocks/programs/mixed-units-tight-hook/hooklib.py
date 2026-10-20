# hooklib.py -- t514-search (wave 16, 13 Oct 2026).  Own code (nothing imported from other agents).
# The hook module M = V (x) Lambda^4 V (dim 560) at block size 3, f = 5, over F_p (p = 1 mod 4, i = sqrt(-1) mod p),
# N = span{e_x (x) e_T : x in T} (dim 280).  rho(g) = g (x) Lambda^4 g.
# Group elements are kept EXACTLY as 8x8 Gaussian-dyadic matrices (int64 re/im arrays times 2^-s) and mapped to F_p
# only for the module.  Arithmetic mod p uses float64 BLAS: p < 2^20 so a 560-term dot product stays below 2^53.
import itertools, sys
import numpy as np

n8 = 8
S4 = list(itertools.combinations(range(n8), 4))
SIDX = {T: i for i, T in enumerate(S4)}
BASIS = [(x, T) for x in range(n8) for T in S4]
NIDX = np.array([i for i, (x, T) in enumerate(BASIS) if x in T])
OIDX = np.array([i for i, (x, T) in enumerate(BASIS) if x not in T])
DM, DN = len(BASIS), len(NIDX)


def isprime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True


def primes_1mod4_below(lim, k):
    out, q = [], lim - 1
    while len(out) < k:
        if q % 4 == 1 and isprime(q): out.append(q)
        q -= 1
    return out


class Field:
    def __init__(self, p):
        self.p = p
        self.I = next(pow(g, (p - 1) // 4, p) for g in range(2, 1000) if pow(g, (p - 1) // 2, p) == p - 1)
        assert self.I * self.I % p == p - 1
        self.inv2 = pow(2, p - 2, p)
        self.pf = float(p)

    def mm(self, A, B):            # A, B float64 arrays with entries in [0, p)
        return np.fmod(A @ B, self.pf)

    def rank(self, A):
        A = np.array(A, dtype=np.int64) % self.p
        p = self.p
        r = 0
        rows, cols = A.shape
        for c in range(cols):
            if r == rows: break
            nz = np.nonzero(A[r:, c])[0]
            if len(nz) == 0: continue
            piv = r + nz[0]
            if piv != r: A[[r, piv]] = A[[piv, r]]
            A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
            col = A[:, c].copy(); col[r] = 0
            nzr = np.nonzero(col)[0]
            if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
            r += 1
        return r

    def rref_basis(self, A):
        """row-reduced basis of the row space (int64)"""
        A = np.array(A, dtype=np.int64) % self.p
        p = self.p
        r = 0
        rows, cols = A.shape
        piv = []
        for c in range(cols):
            if r == rows: break
            nz = np.nonzero(A[r:, c])[0]
            if len(nz) == 0: continue
            pv = r + nz[0]
            if pv != r: A[[r, pv]] = A[[pv, r]]
            A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
            col = A[:, c].copy(); col[r] = 0
            nzr = np.nonzero(col)[0]
            if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
            piv.append(c)
            r += 1
        return A[:r], piv

    def nullspace(self, A):
        """basis of {v : A v = 0} as columns (int64)"""
        A = np.array(A, dtype=np.int64) % self.p
        p = self.p
        R, piv = self.rref_basis(A)
        cols = A.shape[1]
        free = [c for c in range(cols) if c not in set(piv)]
        out = np.zeros((cols, len(free)), dtype=np.int64)
        for k, fc in enumerate(free):
            out[fc, k] = 1
            for i, pc in enumerate(piv):
                out[pc, k] = (-R[i, fc]) % p
        return out


# ---------------------------------------------------------------- exact 8x8 Gaussian-dyadic matrices -------------
class GM:
    """exact matrix (re + i im) / 2^s with int64 entries, s minimal"""
    __slots__ = ("re", "im", "s", "_key")

    def __init__(self, re, im, s):
        re = np.asarray(re, dtype=np.int64); im = np.asarray(im, dtype=np.int64)
        while s > 0 and not ((re | im) & 1).any():
            re = re >> 1; im = im >> 1; s -= 1
        self.re, self.im, self.s = re, im, s
        self._key = None

    def __matmul__(self, o):
        re = self.re @ o.re - self.im @ o.im
        im = self.re @ o.im + self.im @ o.re
        return GM(re, im, self.s + o.s)

    def key(self):
        if self._key is None:
            self._key = self.re.tobytes() + self.im.tobytes() + bytes([self.s])
        return self._key

    def modp(self, F):
        p = F.p
        sc = pow(F.inv2, self.s, p)
        m = (self.re % p + (self.im % p) * F.I) % p
        return (m * sc) % p

    def is_monomial(self):
        nz = (self.re != 0) | (self.im != 0)
        return bool((nz.sum(axis=0) == 1).all() and (nz.sum(axis=1) == 1).all())

    def inv(self):          # only for monomials with unit-modulus entries (phases 1, i, -1, -i) and s = 0
        assert self.s == 0 and self.is_monomial()
        # inverse = conjugate transpose
        return GM(self.re.T.copy(), -self.im.T.copy(), 0)

    def cplx(self):
        return (self.re + 1j * self.im) / (2.0 ** self.s)


def gm_eye():
    return GM(np.eye(8, dtype=np.int64), np.zeros((8, 8), np.int64), 0)


def gm_dir(z):          # dir(z) = (1+i)/2 I + (1-i)/2 X^z
    re = np.zeros((8, 8), np.int64); im = np.zeros((8, 8), np.int64)
    for x in range(8):
        re[x, x] += 1; im[x, x] += 1
        re[x, x ^ z] += 1; im[x, x ^ z] -= 1
    return GM(re, im, 1)


def gm_mono(perm, ph):  # M[perm[x], x] = i^ph[x]   (column x goes to row perm[x])
    re = np.zeros((8, 8), np.int64); im = np.zeros((8, 8), np.int64)
    for x in range(8):
        c = [(1, 0), (0, 1), (-1, 0), (0, -1)][ph[x] % 4]
        re[perm[x], x] = c[0]; im[perm[x], x] = c[1]
    return GM(re, im, 0)


def gm_shift(z):
    return gm_mono([x ^ z for x in range(8)], [0] * 8)


def gm_phase(w, c=0):   # diag( i^{ (w.x + c) mod 2 } )   (the engine's phase letter)
    return gm_mono(list(range(8)), [((bin(w & x).count('1') + c) & 1) for x in range(8)])


def gm_K():
    return gm_dir(1) @ gm_dir(2) @ gm_dir(4)


# ---------------------------------------------------------------- the hook module mod p -------------------------------
_PERM_SIGN = []
for pi in itertools.permutations(range(4)):
    sg = 1
    for a in range(4):
        for b in range(a + 1, 4):
            if pi[a] > pi[b]: sg = -sg
    _PERM_SIGN.append((pi, sg))


def lambda4(gp, F):
    """Lambda^4 of an 8x8 matrix mod p via 4x4 minors (int64 in/out); vectorised over all 70x70 minors"""
    p = F.p
    g = gp.astype(np.int64) % p
    U = np.array(S4)  # 70 x 4
    # sub[i, j, a, b] = g[U[i,a], U[j,b]]
    sub = g[U[:, None, :, None], U[None, :, None, :]]   # 70,70,4,4
    tot = np.zeros((70, 70), dtype=np.int64)
    for pi, sg in _PERM_SIGN:
        t = np.ones((70, 70), dtype=np.int64)
        for a in range(4):
            t = (t * sub[:, :, a, pi[a]]) % p
        tot = (tot + sg * t) % p
    return tot % p


def rho(gp, F):
    """rho(g) on M as float64 mod p (560 x 560)"""
    L = lambda4(gp, F)
    return np.fmod(np.kron(gp.astype(np.float64), L.astype(np.float64)), F.pf)


class Hook:
    def __init__(self, p):
        self.F = F = Field(p)
        self.cache = {}
        Kp = gm_K().modp(F)
        RK = rho(Kp, F)
        # Z = N cap K^-1 N = {n in N : pi_O rho(K) n = 0}
        A = RK[np.ix_(OIDX, NIDX)]
        Zb = F.nullspace(A)                           # DN x 56 (coordinates in N)
        self.Z = Zb.astype(np.float64)
        self.dimZ = Zb.shape[1]
        self.RK = RK

    def rho(self, g):
        k = g.key()
        r = self.cache.get(k)
        if r is None:
            r = rho(g.modp(self.F), self.F)
            if len(self.cache) > 150: self.cache.clear()
            self.cache[k] = r
        return r

    def out_of_N(self, R):          # pi_O R |_N : the 280 x 280 block
        return R[np.ix_(OIDX, NIDX)]
