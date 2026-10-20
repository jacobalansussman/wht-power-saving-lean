"""rt.py -- tb4-nonlinear (wave 16, 13 Oct 2026). Shortest routes to the kernel frame, block size m, in the
WIDENED game (free = every monomial matrix with Z4 entries) and the AFFINE (engine) game.

Facts used [proved in the NOTE]: on a shortest route every row of a depth-d element has support exactly 2^d, and a
backward step from depth d to d-1 is a perfect matching of the rows into pairs (p,q) of equal support whose
entrywise ratio takes two values {v,-v}, v in Z4, each on half the support (b4-search paths.py). A coset B g is
its multiset of rows up to Z4 phase (frame key). Matrices are complex128 with dyadic entries (exact), snapped."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sys, itertools, hashlib
import numpy as np
np.seterr(all="ignore")

A_, B_ = (1 + 1j) / 2, (1 - 1j) / 2
SC = float(2 ** 16)
Z4 = [1, 1j, -1, -1j]


class Game:
    def __init__(self, m):
        self.m = m
        self.n = n = 2 ** m
        self.I = np.eye(n, dtype=complex)
        K = self.I.copy()
        for j in range(m):
            K = self.snap(K @ self.Dm(1 << j))
        self.K = K
        self.memo = {}

    def snap(self, M):
        R = np.round(M.real * SC) / SC
        Im = np.round(M.imag * SC) / SC
        assert np.max(np.abs(R - M.real)) < 1e-9 and np.max(np.abs(Im - M.imag)) < 1e-9, "off grid"
        return R + 1j * Im

    def Xm(self, z):
        M = np.zeros((self.n, self.n), dtype=complex)
        for x in range(self.n):
            M[x ^ z, x] = 1
        return M

    def Dm(self, z):
        return A_ * self.I + B_ * self.Xm(z)

    # ---- frame keys -------------------------------------------------------------------------------
    def ints(self, M):
        return np.round(M.real * SC).astype(np.int64), np.round(M.imag * SC).astype(np.int64)

    def norm_rows(self, M):
        """multiply each row by the power of i making its first nonzero entry p+iq with p>0, q>=0"""
        R, I = self.ints(M)
        n = self.n
        nz = (R != 0) | (I != 0)
        first = nz.argmax(axis=1)
        p = R[np.arange(n), first]; q = I[np.arange(n), first]
        R2, I2 = R.copy(), I.copy()
        for _ in range(4):
            bad = ~((p > 0) & (q >= 0))
            if not bad.any():
                break
            # multiply bad rows by i: (R + iI) i = -I + iR
            R2[bad], I2[bad] = -I2[bad], R2[bad].copy()
            p[bad], q[bad] = -q[bad], p[bad].copy()
        assert ((p > 0) & (q >= 0)).all()
        return R2, I2

    def fkey(self, M):
        R, I = self.norm_rows(M)
        rows = sorted(R[i].tobytes() + I[i].tobytes() for i in range(self.n))
        return hashlib.blake2b(b"".join(rows), digest_size=20).digest()   # 160-bit digest of the frame

    # ---- halving steps ----------------------------------------------------------------------------
    def compat(self, g):
        n = self.n
        nzm = np.abs(g) > 1e-12
        sup = [tuple(np.nonzero(nzm[i])[0]) for i in range(n)]
        E = {}
        for p in range(n):
            for q in range(p + 1, n):
                if sup[p] != sup[q]:
                    continue
                idx = list(sup[p])
                ratio = g[q, idx] / g[p, idx]
                v = ratio[0]
                vz = None
                for z in Z4:
                    if abs(v - z) < 1e-9:
                        vz = z
                if vz is None:
                    continue
                plus = np.abs(ratio - vz) < 1e-9
                minus = np.abs(ratio + vz) < 1e-9
                if not (plus | minus).all() or plus.sum() != minus.sum():
                    continue
                E[(p, q)] = vz
        return E

    def matchings(self, E, limit=None):
        n = self.n
        adj = {p: [] for p in range(n)}
        for (p, q) in E:
            adj[p].append(q); adj[q].append(p)
        out = []
        used = [False] * n

        def rec(cur):
            if limit is not None and len(out) >= limit:
                return
            try:
                p = used.index(False)
            except ValueError:
                out.append(list(cur)); return
            used[p] = True
            for q in adj[p]:
                if not used[q]:
                    used[q] = True
                    cur.append((p, q) if p < q else (q, p))
                    rec(cur)
                    cur.pop()
                    used[q] = False
            used[p] = False
        rec([])
        return out

    def back_step(self, g, M, E):
        h = np.zeros_like(g)
        for k, (p, q) in enumerate(M):
            v = E[(p, q)]
            kq = 1j / v
            rp, rq = g[p], kq * g[q]
            h[2 * k] = B_ * rp + A_ * rq
            h[2 * k + 1] = A_ * rp + B_ * rq
        return self.snap(h)

    def support(self, g):
        s = (np.abs(g) > 1e-12).sum(axis=1)
        return s

    def children(self, g):
        """all halving backward children (as (matching, element)) of g, deduplicated by frame key"""
        E = self.compat(g)
        out = {}
        for M in self.matchings(E):
            h = self.back_step(g, M, E)
            k = self.fkey(h)
            if k not in out:
                out[k] = (M, h)
        return out

    def continues(self, g, lev):
        """can coset Bg (rows of support 2^lev) reach B by lev halving steps?"""
        if lev == 0:
            return True
        k = self.fkey(g)
        if k in self.memo:
            return self.memo[k]
        res = False
        E = self.compat(g)
        for M in self.matchings(E):
            h = self.back_step(g, M, E)
            if self.continues(h, lev - 1):
                res = True
                break
        self.memo[k] = res
        return res

    def nroutes(self, g, lev, memo):
        """number of halving chains (sequences of distinct cosets) from Bg down to B"""
        if lev == 0:
            return 1
        k = self.fkey(g)
        if k in memo:
            return memo[k]
        tot = 0
        for kk, (M, h) in self.children(g).items():
            tot += self.nroutes(h, lev - 1, memo)
        memo[k] = tot
        return tot


def all_matchings(n):
    out = []
    def rec(rem, cur):
        if not rem:
            out.append(list(cur)); return
        p = rem[0]
        for q in rem[1:]:
            cur.append((p, q))
            rec([r for r in rem if r != p and r != q], cur)
            cur.pop()
    rec(list(range(n)), [])
    return out


def gl_elements(m):
    """all invertible m x m bit matrices, as column-lists (images of e_j)"""
    out = []
    n = 2 ** m
    for cols in itertools.product(range(1, n), repeat=m):
        # independence: span size 2^m
        span = {0}
        for c in cols:
            span |= {s ^ c for s in span}
        if len(span) == n:
            out.append(cols)
    return out


def agl_perms(m):
    n = 2 ** m
    G = []
    for cols in gl_elements(m):
        lin = [0] * n
        for x in range(n):
            y = 0
            for j in range(m):
                if (x >> j) & 1:
                    y ^= cols[j]
            lin[x] = y
        for c in range(n):
            G.append([lin[x] ^ c for x in range(n)])
    return np.array(G, dtype=np.int64)
