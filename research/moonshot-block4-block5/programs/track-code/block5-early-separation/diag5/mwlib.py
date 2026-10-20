"""mwlib.py -- mw-big56: a small exact engine for formal words in the ENGINE game at block m (m <= 6).
Group elements are exact 2^m x 2^m matrices on the dyadic grid (entries in Z[i]/2^j), interned by their bytes, so
different scalar multiples are different group elements (the formal group algebra C[Gamma]).
Free letters: affine address maps, Clifford Z4 diagonal phases i^(a.x + 2 x^T B x), Paulis (all monomial, Z4 entries).
Direction letter D = dir(e1) = (1+i)/2 I + (1-i)/2 X_e1 (as in bm.py / fgeo.py); K = dir(e1)...dir(e_m).
Layer of g = log2 max row support (= frame distance from the raw input L0); dK(g) = log2 max row support of K g^-1.
Generators: dict {(l, gid): Fraction}.  Final test: closure over aligned free images (bm.py's method) + exact solve,
with a LAYERED variant: the deepest layer d such that the identity holds on all layers >= d."""
import random
from fractions import Fraction
import numpy as np
np.seterr(all="ignore")
SC = 2.0 ** 14
Z4 = [1, 1j, -1, -1j]

class Eng:
    def __init__(self, m, seed=1):
        self.m = m; self.n = n = 1 << m; self.rng = random.Random(seed)
        self.I = np.eye(n, dtype=complex)
        self.MATS = []; self.IDX = {}; self._lay = {}; self._dk = {}; self._frm = {}; self._mono = {}
        self.Xs = {z: self._X(z) for z in range(n)}
        self.D = self.dirm(1)
        K = self.I
        for j in range(m):
            K = self.snap(K @ self.dirm(1 << j))
        self.K = K; self.KH = np.conj(K).T
        self.ID = self.intern(self.I); self.KID = self.intern(K)

    # ---- exact matrices
    def snap(self, M):
        R = np.round(M.real * SC) / SC; Im = np.round(M.imag * SC) / SC
        assert np.max(np.abs(R - M.real)) < 1e-9 and np.max(np.abs(Im - M.imag)) < 1e-9, "off grid"
        return R + 1j * Im

    def _X(self, z):
        M = np.zeros((self.n, self.n), complex)
        for x in range(self.n):
            M[x ^ z, x] = 1
        return M

    def dirm(self, z):
        return (1 + 1j) / 2 * self.I + (1 - 1j) / 2 * self.Xs[z]

    @staticmethod
    def key(M):
        return np.round(M.real * SC).astype(np.int64).tobytes() + np.round(M.imag * SC).astype(np.int64).tobytes()

    def intern(self, M):
        k = self.key(M); i = self.IDX.get(k)
        if i is None:
            i = len(self.MATS); self.IDX[k] = i; self.MATS.append(M)
        return i

    def mulid(self, M, gid):
        return self.intern(self.snap(M @ self.MATS[gid]))

    # ---- invariants
    def layer(self, gid):
        v = self._lay.get(gid)
        if v is None:
            v = int(round(np.log2((np.abs(self.MATS[gid]) > 1e-9).sum(axis=1).max()))); self._lay[gid] = v
        return v

    def dK(self, gid):
        v = self._dk.get(gid)
        if v is None:
            h = self.snap(self.K @ np.conj(self.MATS[gid]).T)
            v = int(round(np.log2((np.abs(h) > 1e-9).sum(axis=1).max()))); self._dk[gid] = v
        return v

    def is_mono(self, M):
        nz = np.abs(M) > 1e-9
        if not (nz.sum(axis=1) == 1).all() or not (nz.sum(axis=0) == 1).all():
            return False
        return all(min(abs(x - z) for z in Z4) < 1e-9 for x in M[nz])

    def frame(self, gid):
        v = self._frm.get(gid)
        if v is None:
            M = self.MATS[gid]; rows = []
            for i in range(self.n):
                r = M[i]; nz = np.nonzero(np.abs(r) > 1e-9)[0]; c = r[nz[0]]
                for z in Z4:
                    a = np.angle(c * z)
                    if -1e-9 <= a < np.pi / 2 - 1e-9:
                        r = self.snap((r * z)[None, :])[0]; break
                rows.append(self.key(r[None, :]))
            v = b"".join(sorted(rows)); self._frm[gid] = v
        return v

    # ---- free letters (engine)
    def affine(self, A=None, c=None):
        m, n, rng = self.m, self.n, self.rng
        if A is None:
            while True:
                A = [rng.randrange(1, n) for _ in range(m)]
                span = {0}; ok = True
                for a in A:
                    if a in span:
                        ok = False; break
                    span |= {s ^ a for s in span}
                if ok:
                    break
        if c is None:
            c = rng.randrange(n)
        M = np.zeros((n, n), complex)
        for x in range(n):
            y = c
            for j in range(m):
                if (x >> j) & 1:
                    y ^= A[j]
            M[y, x] = 1
        return M

    def phase(self, a=None, B=None):
        m, n, rng = self.m, self.n, self.rng
        if a is None:
            a = [rng.randrange(4) for _ in range(m)]
        if B is None:
            B = {(j, k): rng.randrange(2) for j in range(m) for k in range(j + 1, m)}
        d = []
        for x in range(n):
            b = [(x >> j) & 1 for j in range(m)]
            v = sum(a[j] * b[j] for j in range(m)) + 2 * sum(B.get((j, k), 0) * b[j] * b[k] for j in range(m) for k in range(j + 1, m))
            d.append(Z4[v % 4])
        return np.diag(d)

    def rand_free(self):
        return self.affine() @ self.phase() @ self.Xs[self.rng.randrange(self.n)]

    # ---- generators
    def inp(self, l):
        return {(l, self.ID): Fraction(1)}

    def lmul(self, M, gen):
        out = {}
        for (l, g), c in gen.items():
            t = (l, self.mulid(M, g)); out[t] = out.get(t, 0) + c
        return {t: c for t, c in out.items() if c != 0}

    @staticmethod
    def add(*pairs):
        out = {}
        for c0, gen in pairs:
            for t, c in gen.items():
                out[t] = out.get(t, 0) + c0 * c
        return {t: c for t, c in out.items() if c != 0}

    def move(self, content):
        return self.lmul(self.D, content)

    # ---- final test: closure (aligned single free images of generators) + layered exact solve
    def closure(self, gens, target, cap=4000):
        index = {}
        for i, g in enumerate(gens):
            for (l, gid) in g:
                if self.layer(gid) > 0:
                    index.setdefault((self.frame(gid), l), []).append((i, gid))
        rows = {target: 0}; cols = {}; ent = []; queue = [target]
        while queue:
            (lt, gt) = queue.pop()
            for (i, gterm) in index.get((self.frame(gt), lt), ()):
                bm = self.snap(self.MATS[gt] @ np.conj(self.MATS[gterm]).T)
                if not self.is_mono(bm):
                    continue
                b = self.intern(bm)
                if (i, b) in cols:
                    continue
                c = len(cols); cols[(i, b)] = c
                if c >= cap:
                    return None
                for (l2, g2), co in gens[i].items():
                    g3 = self.mulid(bm, g2)
                    if self.layer(g3) == 0:
                        continue          # absorbed by free images of the inputs
                    t3 = (l2, g3)
                    r3 = rows.get(t3)
                    if r3 is None:
                        r3 = len(rows); rows[t3] = r3; queue.append(t3)
                    ent.append((r3, c, co))
        return rows, cols, ent

    @staticmethod
    def consistent(nr, nc, ent, rowset):
        """exact over Q: is e_0 (restricted to rowset) in the column span (restricted to rowset)?"""
        ridx = {r: k for k, r in enumerate(sorted(rowset))}
        A = [[Fraction(0)] * (nc + 1) for _ in range(len(ridx))]
        for (r, c, co) in ent:
            if r in ridx:
                A[ridx[r]][c] += co
        if 0 in ridx:
            A[ridx[0]][nc] = Fraction(1)
        nrr = len(A); piv = 0
        for c in range(nc):
            pr = next((r for r in range(piv, nrr) if A[r][c]), None)
            if pr is None:
                continue
            A[piv], A[pr] = A[pr], A[piv]
            inv = 1 / A[piv][c]; A[piv] = [x * inv for x in A[piv]]
            for r in range(nrr):
                if r != piv and A[r][c]:
                    f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[piv])]
            piv += 1
            if piv == nrr:
                break
        return not any(A[r][nc] for r in range(piv, nrr))

    def layered_test(self, gens, l, cap=4000):
        """returns (deepest layer d such that the target's identity holds on every layer >= d, or None if even the
        top fails; 1 means a full unit), number of rows, columns"""
        cl = self.closure(gens, (l, self.KID), cap)
        if cl is None:
            return "cap", 0, 0
        rows, cols, ent = cl
        lay = {r: self.layer(t[1]) for t, r in rows.items()}
        best = None
        for d in range(self.m, 0, -1):
            rs = {r for r, L in lay.items() if L >= d}
            if self.consistent(len(rows), len(cols), ent, rs):
                best = d
            else:
                break
        return best, len(rows), len(cols)
