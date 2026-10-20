# vlib.py -- v-ms-round2 (skeptic), 14 Oct 2026.  Own code, nothing imported from ms-round2.
# Block 3 hook M = V (x) Lambda^4 V (dim 560), N = span{e_x (x) e_T : x in T} (280), arithmetic mod p (p = 1 mod 4).
# Lambda^4 is computed by iterated wedges (not by minors), monomials by direct index maps.
import itertools, numpy as np

P_DEFAULT = 786433          # 3 * 2^18 + 1, prime, = 1 mod 4 (not one of ms-round2's primes)


class Fp:
    def __init__(self, p=P_DEFAULT):
        self.p = p
        a = 2
        while pow(a, (p - 1) // 2, p) != p - 1: a += 1
        self.i = pow(a, (p - 1) // 4, p)
        assert self.i * self.i % p == p - 1
        self.h = (p + 1) // 2           # 1/2
        self.pf = float(p)

    def mm(self, A, B):
        return np.fmod(np.asarray(A, np.float64) @ np.asarray(B, np.float64), self.pf)

    def rref(self, A):
        p = self.p
        A = np.array(A, dtype=np.float64) % p
        rows, cols = A.shape
        r = 0; piv = []
        for c in range(cols):
            if r == rows: break
            nz = np.nonzero(A[r:, c])[0]
            if len(nz) == 0: continue
            k = r + nz[0]
            if k != r: A[[r, k]] = A[[k, r]]
            inv = pow(int(A[r, c]), p - 2, p)
            A[r] = np.fmod(A[r] * inv, self.pf)
            col = A[:, c].copy(); col[r] = 0
            nzr = np.nonzero(col)[0]
            if len(nzr):
                A[nzr] = np.fmod(A[nzr] - np.fmod(col[nzr, None] * A[r][None, :], self.pf) + self.pf, self.pf)
            piv.append(c); r += 1
        return A[:r], piv

    def rank(self, A):
        A = np.asarray(A)
        if A.size == 0: return 0
        if A.shape[0] > A.shape[1]: A = A.T
        return len(self.rref(A)[1])

    def null(self, A):
        """columns spanning {v : A v = 0}"""
        A = np.asarray(A)
        cols = A.shape[1]
        R, piv = self.rref(A)
        ps = set(piv)
        free = [c for c in range(cols) if c not in ps]
        out = np.zeros((cols, len(free)))
        for k, fc in enumerate(free):
            out[fc, k] = 1
            for j, pc in enumerate(piv):
                out[pc, k] = (-R[j, fc]) % self.p
        return out


S4 = list(itertools.combinations(range(8), 4))
SI = {T: k for k, T in enumerate(S4)}
NIDX = np.array([x * 70 + k for x in range(8) for k, T in enumerate(S4) if x in T])
OIDX = np.array([x * 70 + k for x in range(8) for k, T in enumerate(S4) if x not in T])


def lam4(F, g):
    """Lambda^4 g (70 x 70) by iterated wedges of columns"""
    p = F.p
    g = np.asarray(g, dtype=np.int64) % p
    L = np.zeros((70, 70))
    for ci, T in enumerate(S4):
        cur = {(): 1}
        for t in T:
            v = g[:, t]
            new = {}
            for S, a in cur.items():
                for j in range(8):
                    if v[j] == 0 or j in S: continue
                    sg = -1 if (sum(1 for s in S if s > j) % 2) else 1
                    S2 = tuple(sorted(S + (j,)))
                    new[S2] = (new.get(S2, 0) + sg * a * int(v[j])) % p
            cur = new
        for S, a in cur.items(): L[SI[S], ci] = a % p
    return L


def rho(F, g):
    return np.fmod(np.kron(np.asarray(g, np.float64) % F.p, lam4(F, g)), F.pf)


def dirz(F, z):
    p = F.p
    a = (1 + F.i) * F.h % p; b = (1 - F.i) * F.h % p
    M = np.zeros((8, 8), dtype=np.int64)
    for x in range(8):
        M[x, x] = (M[x, x] + a) % p
        M[x ^ z, x] = (M[x ^ z, x] + b) % p
    return M


def mono(F, perm, ph):
    """h e_x = i^ph[x] e_perm[x]"""
    M = np.zeros((8, 8), dtype=np.int64)
    for x in range(8): M[perm[x], x] = pow(F.i, ph[x] % 4, F.p)
    return M


def mmod(F, *Ms):
    out = np.asarray(Ms[0], np.int64) % F.p
    for M in Ms[1:]: out = (out @ (np.asarray(M, np.int64) % F.p)) % F.p
    return out


def rho_mono_index(perm, ph):
    """rho(h) for a monomial h on M as (target index array, phase exponent array): rho(h) e_j = i^e[j] e_tgt[j]"""
    tgt = np.zeros(560, dtype=np.int64); ex = np.zeros(560, dtype=np.int64)
    for x in range(8):
        for k, T in enumerate(S4):
            img = [perm[t] for t in T]
            # sign of sorting img
            inv = sum(1 for a in range(4) for b in range(a + 1, 4) if img[a] > img[b])
            e = ph[x] + sum(ph[t] for t in T) + (2 if inv % 2 else 0)
            tgt[x * 70 + k] = perm[x] * 70 + SI[tuple(sorted(img))]
            ex[x * 70 + k] = e % 4
    return tgt, ex


def apply_mono(F, tgt, ex, X):
    """rho(h) X for X with 560 rows"""
    out = np.zeros_like(X)
    ipow = np.array([pow(F.i, k, F.p) for k in range(4)], dtype=np.float64)
    out[tgt] = np.fmod(X * ipow[ex][:, None], F.pf)
    return out
