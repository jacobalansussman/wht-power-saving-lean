"""bm.py -- tb4-beam (wave 16, 13 Oct 2026). Guided exact search in the WIDENED game (free = every monomial matrix with
entries in {1,i,-1,-i}) at block size m = 3 or 4.

Game (GAME.md, register form): generators start as the inputs x_l; a move makes u = D.c with D = dir(e1) (rank 1) or
dir(e1)dir(e2) (a rank-2 block), c = a Q[B']-combination of earlier generators (left free factors are absorbed by later
contents); end: [K]x_l in the Q[B']-span of the generators, every l.  Exact part: formal sums of 2^m x 2^m dyadic
matrices (complex128, exact, snapped and asserted), final test = closure + exact solve over Q (s-alg / b4-search eng4),
free last move test = f-reduce pfinal / b4-search free4 (ported to any m, and allowing every term with K g^-1 one move
from a monomial).  Score part: Theorem R potential in the hook module V (x) Lambda^k V, N = {e_a (x) e_S : a in S},
computed mod p = 1048573 (float64 BLAS, ranks in C), restricted to N (the contents' images lie in N on S_t).
INSIDE gain of a move = drop of dim(S_t + E^W), E = N cap K^-1 N; needed in total W*dem; waste = cap - inside gain."""
import sys, os, itertools, random, ctypes, time
PYLIB = ''
if os.path.isdir(PYLIB):
    sys.path.insert(0, PYLIB)
import numpy as np
from fractions import Fraction
np.seterr(all="ignore")
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- mod p linear algebra
P = 1048573
IP = next(pow(g, (P - 1) // 4, P) for g in range(2, 500) if pow(g, (P - 1) // 2, P) == P - 1)
assert IP * IP % P == P - 1
_lib = ctypes.CDLL(os.path.join(HERE, 'libmp_mac.so' if sys.platform == 'darwin' else 'libmp_lin.so'))
_lib.rref.restype = ctypes.c_long
_lib.rref.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_uint64, ctypes.c_void_p, ctypes.c_long]

def mmp(A, B):
    """A @ B mod P, float64 exact while inner dim < 8192"""
    assert A.shape[1] < 8192
    return np.mod(np.asarray(A, dtype=np.float64) @ np.asarray(B, dtype=np.float64), P)

def rref(A, full=True):
    A = np.ascontiguousarray(np.mod(A, P).astype(np.uint64))
    R, C = A.shape
    piv = np.zeros(max(1, min(R, C)), dtype=np.int64)
    r = _lib.rref(A.ctypes.data, R, C, P, piv.ctypes.data, 1 if full else 0)
    return A[:r].astype(np.float64), list(piv[:r])

RNG = np.random.default_rng(12345)
def rank(A):
    if A.size == 0:
        return 0
    R, C = A.shape
    if C > R + 16:                       # compress columns (rank <= R)
        A = mmp(A, RNG.integers(0, P, size=(C, R + 8)).astype(np.float64))
    return len(rref(A, full=False)[1])

def nullspace(A):
    """columns basis of {x : A x = 0}"""
    R, C = A.shape
    Rr, piv = rref(A, full=True)
    free = [c for c in range(C) if c not in set(piv)]
    Nm = np.zeros((C, len(free)))
    for t, fc in enumerate(free):
        Nm[fc, t] = 1
        for i, pc in enumerate(piv):
            Nm[pc, t] = (-Rr[i, fc]) % P
    return Nm

def cplx_to_p(M, SC=4096):
    re = np.round(M.real * SC).astype(np.int64); im = np.round(M.imag * SC).astype(np.int64)
    assert np.allclose(re / SC, M.real) and np.allclose(im / SC, M.imag)
    inv = pow(SC, P - 2, P)
    return np.mod((re % P + (im % P) * IP % P) % P * inv, P).astype(np.int64)

# ---------------------------------------------------------------- exact matrices / group ring
A_, B_ = (1 + 1j) / 2, (1 - 1j) / 2
Z4 = [1, 1j, -1, -1j]
SC = 4096.0

class Game:
    def __init__(self, m):
        self.m = m; self.n = n = 1 << m
        self.I = np.eye(n, dtype=complex)
        self.MATS = []; self.IDX = {}; self.FREE = {}; self.FRM = {}
        self.D = {z: self.Dm(z) for z in range(1, n)}
        K = self.I
        for j in range(m):
            K = self.mul(K, self.D[1 << j])
        self.K = K
        self.KID = self.intern(K); self.ID = self.intern(self.I)
        self.D1 = self.D[1]; self.D12 = self.mul(self.D[1], self.D[2])

    def snap(self, M):
        R = np.round(M.real * SC) / SC; Im = np.round(M.imag * SC) / SC
        assert np.max(np.abs(R - M.real)) < 1e-9 and np.max(np.abs(Im - M.imag)) < 1e-9, "not on grid"
        return R + 1j * Im

    def mul(self, *Ms):
        out = Ms[0]
        for M in Ms[1:]:
            out = self.snap(out @ M)
        return out

    def Xm(self, z):
        M = np.zeros((self.n, self.n), dtype=complex)
        for x in range(self.n):
            M[x ^ z, x] = 1
        return M

    def Dm(self, z):
        return A_ * self.I + B_ * self.Xm(z)

    @staticmethod
    def key(M):
        return (np.round(M.real * SC).astype(np.int64).tobytes() + np.round(M.imag * SC).astype(np.int64).tobytes())

    def intern(self, M):
        k = self.key(M)
        i = self.IDX.get(k)
        if i is None:
            i = len(self.MATS); self.IDX[k] = i; self.MATS.append(M)
        return i

    @staticmethod
    def is_monomial(M):
        nz = np.abs(M) > 1e-9
        if not (nz.sum(axis=1) == 1).all() or not (nz.sum(axis=0) == 1).all():
            return False
        v = M[nz]
        return all(min(abs(x - z) for z in Z4) < 1e-9 for x in v)

    def is_free(self, gid):
        v = self.FREE.get(gid)
        if v is None:
            v = self.is_monomial(self.MATS[gid]); self.FREE[gid] = v
        return v

    def row_phase_norm(self, r):
        nz = np.nonzero(np.abs(r) > 1e-9)[0]
        v = r[nz[0]]
        for z in Z4:
            w = v * z
            ang = np.angle(w)
            if -1e-9 <= ang < np.pi / 2 - 1e-9:
                return self.snap((r * z)[None, :])[0]
        raise AssertionError

    def frame_key(self, M):
        rows = [self.row_phase_norm(M[i]) for i in range(self.n)]
        ks = sorted((np.round(r.real * SC).astype(np.int64).tobytes() + np.round(r.imag * SC).astype(np.int64).tobytes()) for r in rows)
        return b"".join(ks)

    def frm(self, gid):
        v = self.FRM.get(gid)
        if v is None:
            v = self.frame_key(self.MATS[gid]); self.FRM[gid] = v
        return v

    @staticmethod
    def rowsupp(M):
        return (np.abs(M) > 1e-9).sum(axis=1)

    def dist(self, gid):
        """lower bound on the number of moves from frame B'g to the kernel frame: log2 row support of K g^-1"""
        h = self.mul(self.K, np.conj(self.MATS[gid]).T)
        return int(round(np.log2(self.rowsupp(h).max())))

    def mono(self, sigma, q):
        """monomial e_y -> i^q[y] e_sigma[y]"""
        M = np.zeros((self.n, self.n), dtype=complex)
        for y in range(self.n):
            M[sigma[y], y] = Z4[q[y] % 4]
        return M

    @staticmethod
    def mono_of(M):
        n = M.shape[0]
        sigma = [0] * n; q = [0] * n
        for y in range(n):
            x = int(np.argmax(np.abs(M[:, y])))
            sigma[y] = x
            q[y] = min(range(4), key=lambda t: abs(M[x, y] - Z4[t]))
        return sigma, q

    # ---- generators: dict {(l, gid): Fraction}
    def inp(self, l):
        return {(l, self.ID): Fraction(1)}

    def lmul(self, M, gen):
        out = {}
        for (l, g), c in gen.items():
            t = (l, self.intern(self.mul(M, self.MATS[g])))
            out[t] = out.get(t, 0) + c
        return {t: c for t, c in out.items() if c != 0}

    @staticmethod
    def add(*pairs):
        out = {}
        for (c0, gen) in pairs:
            for t, c in gen.items():
                out[t] = out.get(t, 0) + c0 * c
        return {t: c for t, c in out.items() if c != 0}

    # ---- final test (b4-search eng4, generic m)
    def closure(self, gens, target, cap, E=None):
        index = {}
        for i, g in enumerate(gens):
            for (ll, gid) in g:
                index.setdefault((self.frm(gid), ll), []).append((i, gid))
        rows = {target: 0}; cols = {}; ent = []; queue = [target]
        Einv = None if E is None else np.conj(E).T
        def freeish(g3):
            if self.is_free(g3):
                return True
            return Einv is not None and self.is_monomial(self.mul(Einv, self.MATS[g3]))
        while queue:
            (lt, gt) = queue.pop()
            cands = []
            for (i, gterm) in index.get((self.frm(gt), lt), ()):
                cands.append((0, i, self.mul(self.MATS[gt], np.conj(self.MATS[gterm]).T)))
            if E is not None:
                fk = self.frame_key(self.mul(Einv, self.MATS[gt]))
                for (i, gterm) in index.get((fk, lt), ()):
                    cands.append((1, i, self.mul(Einv, self.MATS[gt], np.conj(self.MATS[gterm]).T)))
            for (kind, i, bm) in cands:
                if not self.is_monomial(bm):
                    continue
                b = self.intern(bm)
                ck = (kind, i, b)
                if ck in cols:
                    continue
                c = len(cols); cols[ck] = c
                if c >= cap:
                    return None
                pre = self.mul(E, bm) if kind == 1 else bm
                for (l2, g2), co in gens[i].items():
                    g3 = self.intern(self.mul(pre, self.MATS[g2]))
                    if freeish(g3):
                        continue
                    t3 = (l2, g3)
                    r3 = rows.get(t3)
                    if r3 is None:
                        r3 = len(rows); rows[t3] = r3; queue.append(t3)
                    ent.append((r3, c, co))
        return rows, cols, ent

    @staticmethod
    def solve_modp(nr, nc, ent):
        if nc == 0:
            return False
        A = np.zeros((nr, nc + 1))
        for (r, c, co) in ent:
            A[r, c] = (A[r, c] + (co.numerator % P) * pow(co.denominator % P, P - 2, P)) % P
        A[0, nc] = 1
        Rr, piv = rref(A, full=False)
        return nc not in piv

    @staticmethod
    def solve_Q(nr, nc, ent):
        A = [[Fraction(0)] * (nc + 1) for _ in range(nr)]
        for (r, c, co) in ent:
            A[r][c] += co
        A[0][nc] = Fraction(1)
        piv, pcols = 0, []
        for c in range(nc):
            pr = next((r for r in range(piv, nr) if A[r][c]), None)
            if pr is None:
                continue
            A[piv], A[pr] = A[pr], A[piv]
            inv = 1 / A[piv][c]
            A[piv] = [x * inv for x in A[piv]]
            for r in range(nr):
                if r != piv and A[r][c]:
                    f = A[r][c]
                    A[r] = [x - f * y for x, y in zip(A[r], A[piv])]
            pcols.append(c); piv += 1
            if piv == nr:
                break
        if any(A[r][nc] for r in range(piv, nr)):
            return None
        x = [Fraction(0)] * nc
        for r, c in enumerate(pcols):
            x[c] = A[r][nc]
        return x

    def span_test(self, gens, l, cap=600, E=None):
        cl = self.closure(gens, (l, self.KID), cap, E)
        if cl is None:
            return "cap", None
        rows, cols, ent = cl
        if not self.solve_modp(len(rows), len(cols), ent):
            return "no", None
        x = self.solve_Q(len(rows), len(cols), ent)
        if x is None:
            return "no", None
        inv = {c: k for k, c in cols.items()}
        return "yes", [(inv[c], x[c]) for c in range(len(x)) if x[c]]

    # ---- halving steps (b4-search paths.py, generic n)
    def compat(self, g):
        n = self.n
        sup = [tuple(np.nonzero(np.abs(g[i]) > 1e-9)[0]) for i in range(n)]
        E = {}
        for p_ in range(n):
            for q_ in range(p_ + 1, n):
                if sup[p_] != sup[q_]:
                    continue
                ratio = g[q_, list(sup[p_])] / g[p_, list(sup[p_])]
                v = ratio[0]
                if min(abs(v - z) for z in Z4) > 1e-9:
                    continue
                plus = np.abs(ratio - v) < 1e-9; minus = np.abs(ratio + v) < 1e-9
                if not (plus | minus).all() or plus.sum() != minus.sum():
                    continue
                E[(p_, q_)] = v
        return E

    def rand_matchings(self, E, rng, want, tries):
        """distinct random perfect matchings of the compat graph (randomized DFS)"""
        n = self.n
        adj = {p_: [] for p_ in range(n)}
        for (a, b) in E:
            adj[a].append(b); adj[b].append(a)
        out = {}
        for _ in range(tries):
            if len(out) >= want:
                break
            res = []
            def rec(rem):
                if not rem:
                    return True
                p_ = min(rem, key=lambda x: len([y for y in adj[x] if y in rem]))
                nb = [y for y in adj[p_] if y in rem]
                rng.shuffle(nb)
                for q_ in nb:
                    res.append((min(p_, q_), max(p_, q_)))
                    if rec(rem - {p_, q_}):
                        return True
                    res.pop()
                return False
            if rec(frozenset(range(n))):
                out[tuple(sorted(res))] = 1
        return [list(M) for M in out]

    def back_step(self, g, M, E):
        h = np.zeros_like(g)
        for k, (p_, q_) in enumerate(M):
            v = E[(p_, q_)]
            kq = 1j / v
            rp, rq = g[p_], kq * g[q_]
            h[2 * k] = B_ * rp + A_ * rq
            h[2 * k + 1] = A_ * rp + B_ * rq
        return self.snap(h)

    def one_move(self, h):
        if not (self.rowsupp(h) == 2).all():
            return False
        E = self.compat(h)
        rng = random.Random(1)
        for M in self.rand_matchings(E, rng, 30, 200):
            if self.is_monomial(self.back_step(h, M, E)):
                return True
        return False

    def E_candidates(self, gens, l):
        out = {}
        for g in gens:
            for (ll, gid) in g:
                if ll != l:
                    continue
                h = self.mul(self.K, np.conj(self.MATS[gid]).T)
                if (self.rowsupp(h) == 2).all() and self.one_move(h):
                    out.setdefault(self.frame_key(h.T.copy()), h)
        return list(out.values())

    def free_last_test(self, gens, l, cap=800):
        r, sol = self.span_test(gens, l, cap)
        if r == "yes":
            return "yes0", sol
        Es = self.E_candidates(gens, l)
        if not Es:
            return "none", None
        res = "no"
        for E in Es:
            rr, sol = self.span_test(gens, l, cap, E)
            if rr == "yes":
                return "yes", (E, sol)
            if rr == "cap":
                res = "cap"
        return res, None

    # ---- progress letters: b such that D1 b g is one step nearer the kernel frame
    def continues(self, h, rng, budget=6):
        """does a chain of rank-1 right halvings h -> h b^-1 D1^-1 -> ... reach a monomial? (randomised DFS, memo)"""
        if not hasattr(self, "_cont"):
            self._cont = {}
        k = self.key(h)
        v = self._cont.get(k)
        if v is not None:
            return v
        rs = self.rowsupp(h)
        d = int(round(np.log2(rs.max())))
        if d == 0:
            v = self.is_monomial(h)
        elif not (rs == (1 << d)).all():
            v = False
        else:
            v = False
            E = self.compat(h.T.copy())
            D1i = np.conj(self.D1).T
            for M in self.rand_matchings(E, rng, budget, budget * 3):
                Pi = np.zeros((self.n, self.n), dtype=complex)
                for kk, (p_, q_) in enumerate(M):
                    Pi[2 * kk, p_] = 1; Pi[2 * kk + 1, q_] = 1j / E[(p_, q_)]
                b = np.conj(Pi)
                hn = self.mul(h, np.conj(b).T, D1i)
                if (self.rowsupp(hn) == (1 << (d - 1))).all() and self.continues(hn, rng, budget):
                    v = True
                    break
        self._cont[k] = v
        return v

    def progress_letters(self, gid, rng, want=6, Dmat=None, steps=1):
        """monomials b with dist(D b g) = dist(g) - steps, D = D1 (steps 1) or D1D2 (steps 2). Candidates: the m!
        coordinate permutations (engine moves along coordinate directions) and b from random halving matchings of the
        columns of h = K g^-1 (b = conj(Pi), Pi moving matched columns to (2k, 2k+1)); one-level lookahead (the new
        h must still be halvable unless free). Returns (all valid coordinate ones) + up to `want` random ones."""
        g = self.MATS[gid]
        h = self.mul(self.K, np.conj(g).T)
        d0 = int(round(np.log2(self.rowsupp(h).max())))
        if d0 < steps or not (self.rowsupp(h) == (1 << d0)).all():
            return []
        Dm_ = self.D1 if Dmat is None else Dmat
        def ok(b):
            newg = self.mul(Dm_, b, g)
            hn = self.mul(self.K, np.conj(newg).T)
            dn = int(round(np.log2(self.rowsupp(hn).max())))
            if dn != d0 - steps:
                return False
            return self.continues(hn, rng)       # a full halving chain to a monomial must exist (no dead ends)
        
        out = {}
        frames = set()
        for perm in itertools.permutations(range(self.m)):
            sig = [sum(((x >> j) & 1) << perm[j] for j in range(self.m)) for x in range(self.n)]
            b = self.mono(sig, [0] * self.n)
            if ok(b):
                fk = self.frame_key(self.mul(Dm_, b, g))      # same frame = same move up to a free left factor
                if fk in frames:
                    continue
                frames.add(fk)
                out[self.key(b)] = b
        ncoord = len(out)
        ht = h.T.copy()
        E = self.compat(ht)
        tries = want * 4 if steps == 1 else want * 20
        Ms = self.rand_matchings(E, rng, tries, tries * 4)
        for M in Ms:
            order = list(range(len(M)))
            if steps == 2:
                rng.shuffle(order)
            Pi = np.zeros((self.n, self.n), dtype=complex)
            for k, idx in enumerate(order):
                (p_, q_) = M[idx]
                v = E[(p_, q_)]
                Pi[2 * k, p_] = 1; Pi[2 * k + 1, q_] = 1j / v
            b = np.conj(Pi)
            kb = self.key(b)
            if kb in out:
                continue
            if ok(b):
                fk = self.frame_key(self.mul(Dm_, b, g))
                if fk in frames:
                    continue
                frames.add(fk)
                out[kb] = b
            if len(out) >= ncoord + want:
                break
        return list(out.values())

# ---------------------------------------------------------------- hook module restricted to N
def det_batch(M):
    """batched det mod P of (..., k, k) int64 arrays (Laplace)"""
    k = M.shape[-1]
    if k == 1:
        return M[..., 0, 0] % P
    if k == 2:
        return (M[..., 0, 0] * M[..., 1, 1] - M[..., 0, 1] * M[..., 1, 0]) % P
    out = np.zeros(M.shape[:-2], dtype=np.int64)
    for j in range(k):
        minor = np.delete(np.delete(M, 0, axis=-2), j, axis=-1)
        t = (M[..., 0, j] % P) * det_batch(minor) % P
        out = (out + (t if j % 2 == 0 else -t)) % P
    return out

class HookN:
    def __init__(self, G, k):
        self.G = G; self.k = k; n = G.n
        self.subs = list(itertools.combinations(range(n), k))
        self.sidx = {S: i for i, S in enumerate(self.subs)}
        self.Nb = [(a, S) for S in self.subs for a in S]
        self.Nidx = {b: i for i, b in enumerate(self.Nb)}
        self.dN = len(self.Nb)
        self.comp = [(a, S) for S in self.subs for a in range(n) if a not in S]

    def wedge(self, gp):
        subs = np.array(self.subs)
        sub = gp[subs[:, None, :, None], subs[None, :, None, :]]   # [S', S, i, j] = g[S'_i, S_j]
        return det_batch(sub)

    def rho_cols(self, Mexact):
        """rho(g) restricted to columns N: (rows N, rows comp) as float64 mod P"""
        gp = cplx_to_p(Mexact)
        Wd = self.wedge(gp)
        ncol = self.dN
        ca = np.array([a for (a, S) in self.Nb]); cs = np.array([self.sidx[S] for (a, S) in self.Nb])
        ra = np.array([a for (a, S) in self.Nb]); rs = np.array([self.sidx[S] for (a, S) in self.Nb])
        RN = (gp[ra[:, None], ca[None, :]] * Wd[rs[:, None], cs[None, :]]) % P
        qa = np.array([a for (a, S) in self.comp]); qs = np.array([self.sidx[S] for (a, S) in self.comp])
        RC = (gp[qa[:, None], ca[None, :]] * Wd[qs[:, None], cs[None, :]]) % P
        return RN.astype(np.float64), RC.astype(np.float64)

    def mono_N(self, sigma, q):
        """rho(b) on N for b e_y = i^q[y] e_sigma[y]: returns (dest index array, coefficient array)"""
        dest = np.zeros(self.dN, dtype=np.int64); coef = np.zeros(self.dN)
        for i, (a, S) in enumerate(self.Nb):
            img = [sigma[s] for s in S]
            perm = sorted(range(len(img)), key=lambda t: img[t])
            # sign of the sorting permutation
            sgn = 1; seen = [False] * len(perm)
            for s0 in range(len(perm)):
                if not seen[s0]:
                    L = 0; x = s0
                    while not seen[x]:
                        seen[x] = True; x = perm[x]; L += 1
                    if L % 2 == 0:
                        sgn = -sgn
            ph = (q[a] + sum(q[s] for s in S)) % 4
            c = pow(IP, ph, P) * (1 if sgn > 0 else P - 1) % P
            dest[i] = self.Nidx[(sigma[a], tuple(sorted(img)))]
            coef[i] = c
        return dest, coef

def apply_mono(mn, Z):
    dest, coef = mn
    out = np.empty_like(Z)
    out[dest] = np.mod(coef[:, None] * Z, P)
    return out

class Pot:
    """Theorem-R potential for W arrays, hook module restricted to N, mod P."""
    def __init__(self, G, k, W, letters=("D1",)):
        self.G = G; self.H = HookN(G, k); self.W = W
        dN = self.H.dN
        self.F = {}; self.RNN = {}; self.cap = {}
        mats = {"D1": G.D1, "D12": G.D12}
        for nm in letters:
            RN, RC = self.H.rho_cols(mats[nm])
            cap = rank(RC.T.copy()) if RC.shape[0] > RC.shape[1] else rank(RC)
            R1 = RNG.integers(0, P, size=(cap + 24, RC.shape[0])).astype(np.float64)
            self.F[nm] = mmp(R1, RC) if RC.shape[0] < 8192 else self._bigproj(R1, RC)
            self.RNN[nm] = RN; self.cap[nm] = cap
        RN, RC = self.H.rho_cols(G.K)
        Rk = RNG.integers(0, P, size=(dN + 8, RC.shape[0])).astype(np.float64)
        KK = mmp(Rk, RC) if RC.shape[0] < 8192 else self._bigproj(Rk, RC)
        self.Ebasis = nullspace(KK)                  # dN x e
        self.e = self.Ebasis.shape[1]
        self.dem = dN - self.e

    @staticmethod
    def _bigproj(R1, RC):
        out = np.zeros((R1.shape[0], RC.shape[1]))
        for s in range(0, RC.shape[0], 4096):
            out = np.mod(out + mmp(R1[:, s:s + 4096], RC[s:s + 4096]), P)
        return out

    def root(self):
        W, dN, e = self.W, self.H.dN, self.e
        B = np.eye(W * dN, dtype=np.float32)
        T = np.zeros((W * dN, W * e), dtype=np.float32)
        for l in range(W):
            T[l * dN:(l + 1) * dN, l * e:(l + 1) * e] = self.Ebasis
        return {"B": B, "T": T, "Z": [], "Y": [], "phi": 0}

    def gen_mats(self, st, j):
        """(Z_j, Y_j) for generator j (inputs first)"""
        W, dN = self.W, self.H.dN
        if j < W:
            return st["B"][j * dN:(j + 1) * dN], st["T"][j * dN:(j + 1) * dN]
        return st["Z"][j - W], st["Y"][j - W]

    def content(self, st, terms, onlyY=False):
        """terms: list of (sign, mono_N or None, j)"""
        Zc = None; Yc = None
        for (sg, mn, j) in terms:
            Z, Y = self.gen_mats(st, j)
            if onlyY:
                Z = None
            if mn is not None:
                Y = apply_mono(mn, Y)
                if Z is not None: Z = apply_mono(mn, Z)
            if sg < 0:
                Y = np.mod(-Y, P)
                if Z is not None: Z = np.mod(-Z, P)
            if Z is not None:
                Zc = Z if Zc is None else np.mod(Zc + Z, P)
            Yc = Y if Yc is None else np.mod(Yc + Y, P)
        return Zc, Yc

    def score(self, st, terms, letter="D1"):
        """(raw gain, inside gain). Uses per-node column projections Z_j R (R random s x q, q = rows of F >= cap + 24):
        rank(F Zc) = rank(F Zc R) with high probability since rank(F Zc) <= cap < q."""
        F = self.F[letter]
        if "ZR" not in st:
            s_ = st["B"].shape[1]
            q = max(f.shape[0] for f in self.F.values()) + 8
            if s_ > q:
                R = RNG.integers(0, P, size=(s_, q)).astype(np.float64)
                ng = self.W + len(st["Z"])
                st["ZR"] = [mmp(self.gen_mats(st, j)[0], R) for j in range(ng)]
            else:
                st["ZR"] = [self.gen_mats(st, j)[0] for j in range(self.W + len(st["Z"]))]
        Zc = None
        for (sg, mn, j) in terms:
            Z = st["ZR"][j]
            if mn is not None:
                Z = apply_mono(mn, Z)
            if sg < 0:
                Z = np.mod(-Z, P)
            Zc = Z if Zc is None else np.mod(Zc + Z, P)
        _, Yc = self.content(st, terms, onlyY=True)
        g = rank(mmp(F, Zc)) if Zc.shape[1] else 0
        d = rank(mmp(F, Yc)) if Yc.shape[1] else 0
        return g, g - d

    def commit(self, st, terms, letter="D1"):
        Zc, Yc = self.content(st, terms)
        F = self.F[letter]
        GZ = mmp(F, Zc); GY = mmp(F, Yc)
        Nz = nullspace(GZ) if Zc.shape[1] else np.zeros((0, 0))
        Ny = nullspace(GY) if Yc.shape[1] else np.zeros((0, 0))
        g = Zc.shape[1] - Nz.shape[1]; d = Yc.shape[1] - Ny.shape[1]
        f32 = lambda A: A.astype(np.float32)          # entries < 2^20: exact in float32
        new = {"B": f32(mmp(st["B"], Nz)), "T": f32(mmp(st["T"], Ny)),
               "Z": [f32(mmp(Z, Nz)) for Z in st["Z"]] + [f32(mmp(mmp(self.RNN[letter], Zc), Nz))],
               "Y": [f32(mmp(Y, Ny)) for Y in st["Y"]] + [f32(mmp(mmp(self.RNN[letter], Yc), Ny))],
               "phi": st["phi"] + g - d}
        return new, g, g - d
