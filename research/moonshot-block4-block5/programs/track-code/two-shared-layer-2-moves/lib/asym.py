"""asym.py -- w3-asym: asymmetric programs at block 4 (register form of GAME.md s5), judged by EXACT-structure
representation tests instead of frame sums.
Program: W inputs x_l; move t: u_t = D1 . sum_k alpha_{t,k} b_{t,k} u_{s_{t,k}} (b free letter, s an earlier generator or input).
Coefficients alpha are complex parameters (first summand of each move fixed to 1: projective).
Tests (all NECESSARY for a unit, engine and widened game; W outputs):
  Df(f): diagonal-pattern block of V^(x)f: rows M[x,:]^(x)f of every term; the coordinate projection onto the diagonal
         pattern commutes with every b^(x)f (b monomial), and the image of C[B'] on it lies in M16, so [K]x_l in the span
         forces rows K[x,:]^(x)f (slot l) in the row span of the generators' rows.  f = 2, 3 explicit vectors.
  full f = 2: rep2.Rep2 (D16 + S120 + A120).
Residual = squared distance of the target rows from the row span / target norm, summed over outputs."""
import sys, os, random, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from bm import Game
from rep2 import Rep2

G = Game(4)
n = 16

def bits(x):
    return [(x >> j) & 1 for j in range(4)]

def rand_engine(rng):
    """engine free letter: affine address map (GL(4,2) + shift) times diagonal phase i^(a.x + 2 sum C_jk x_j x_k)"""
    while True:
        A = [rng.randrange(1, n) for _ in range(4)]
        span = {0}; ok = True
        for a in A:
            if a in span: ok = False; break
            span |= {x ^ a for x in span}
        if ok: break
    c = rng.randrange(n)
    sig = []
    for x in range(n):
        y = c
        for j in range(4):
            if (x >> j) & 1: y ^= A[j]
        sig.append(y)
    a = [rng.randrange(4) for _ in range(4)]
    C = {(j, k): rng.randrange(2) for j in range(4) for k in range(j + 1, 4)}
    q = []
    for x in range(n):
        b = bits(x)
        v = sum(a[j] * b[j] for j in range(4)) + 2 * sum(C[(j, k)] * b[j] * b[k] for j in range(4) for k in range(j + 1, 4))
        q.append(v % 4)
    return G.mono(sig, q)

def rand_wide(rng):
    sig = list(range(n)); rng.shuffle(sig)
    return G.mono(sig, [rng.randrange(4) for _ in range(n)])

def dist(gid):
    return G.dist(gid)

class Prog:
    def __init__(self, W):
        self.W = W
        self.terms = [[(l, G.ID)] for l in range(W)]     # per generator: list of (l, gid)
        self.maps = [None] * W                            # per move: list of (param or None, src, si, ri)
        self.np = 0
        self.moves = []                                   # list of summand lists (src, bmatrix, pidx)
    def add_move(self, summ, D=None):
        """summ: list of (src, b) ; first summand coefficient fixed 1, others new params"""
        D = G.D1 if D is None else D
        terms = {}; order = []; maps = []; rec = []
        for k, (src, b) in enumerate(summ):
            pidx = None if k == 0 else self.np + k - 1
            si, ri = [], []
            for j, (l, g) in enumerate(self.terms[src]):
                g3 = G.intern(G.mul(D, b, G.MATS[g]))
                key = (l, g3)
                if key not in terms:
                    terms[key] = len(order); order.append(key)
                si.append(j); ri.append(terms[key])
            maps.append((pidx, src, np.array(si), np.array(ri)))
            rec.append((src, b, pidx))
        self.np += len(summ) - 1
        self.terms.append(order); self.maps.append(maps); self.moves.append(rec)
        return len(self.terms) - 1
    def coefs(self, par):
        out = []
        for i, T in enumerate(self.terms):
            if self.maps[i] is None:
                out.append(np.ones(1, dtype=complex)); continue
            v = np.zeros(len(T), dtype=complex)
            for (pidx, src, si, ri) in self.maps[i]:
                c = 1.0 if pidx is None else par[pidx]
                np.add.at(v, ri, c * out[src][si])
            out.append(v)
        return out
    def describe(self):
        s = []
        for i, T in enumerate(self.terms):
            ds = {}
            for (l, g) in T:
                ds.setdefault(l, []).append(4 - dist(g))
            s.append(f"u{i}:" + ",".join(f"x{l}@{sorted(v)}" for l, v in sorted(ds.items())))
        return " ".join(s)

class DTest:
    """diagonal-pattern tests at f = 2, 3 (explicit), precomputed per term"""
    def __init__(self, prog, fs=(2, 3)):
        self.p = prog; self.fs = fs; self.W = prog.W
        self.rows = {}
        Kr = G.K
        self.T = {f: np.array([self.tens(Kr[x], f) for x in range(n)]) for f in fs}
        for i, T in enumerate(prog.terms):
            for (l, g) in T:
                if g not in self.rows:
                    M = G.MATS[g]
                    self.rows[g] = {f: np.array([self.tens(M[x], f) for x in range(n)]) for f in fs}
    @staticmethod
    def tens(v, f):
        out = v
        for _ in range(f - 1):
            out = np.kron(out, v)
        return out
    def residual(self, par, outputs=None, detail=False):
        cv = self.p.coefs(par)
        tot = 0; det = {}
        outputs = range(self.W) if outputs is None else outputs
        for f in self.fs:
            dim = 16 ** f
            blocks = []
            for i, T in enumerate(self.p.terms):
                X = np.zeros((n, dim * self.W), dtype=complex)
                for j, (l, g) in enumerate(T):
                    X[:, dim * l:dim * (l + 1)] += cv[i][j] * self.rows[g][f]
                blocks.append(X)
            S = np.vstack(blocks)
            U, s, Vh = np.linalg.svd(S, full_matrices=False)
            r = int((s > 1e-9 * s[0]).sum())
            Vr = Vh[:r]
            for l0 in outputs:
                Tm = np.zeros((n, dim * self.W), dtype=complex); Tm[:, dim * l0:dim * (l0 + 1)] = self.T[f]
                R = Tm - (Tm @ Vr.conj().T) @ Vr
                res = float(np.sum(np.abs(R) ** 2) / np.sum(np.abs(Tm) ** 2))
                tot += res; det[(f, l0)] = round(res, 6)
        return (tot, det) if detail else tot

def full2(prog, par):
    R = Rep2(prog.W, G.K)
    cv = prog.coefs(par)
    tl = [[(l, G.MATS[g], cv[i][j]) for j, (l, g) in enumerate(T)] for i, T in enumerate(prog.terms)]
    return [R.residual(tl, l0, detail=True) for l0 in range(prog.W)]

def plain(W):
    P = Prog(W)
    for l in range(W):
        cur = l
        for j in range(4):
            A = np.eye(n, dtype=complex)
            # address map sending e1 -> e_(j+1) direction: D_{e_j} = A D1 A^-1; content letter b = A^-1 applied before
            perm = list(range(n))
            # swap bit 0 and bit j
            def sw(x):
                b0 = x & 1; bj = (x >> j) & 1
                x2 = x & ~1 & ~(1 << j)
                return x2 | (bj) | (b0 << j)
            sig = [sw(x) for x in range(n)]
            Am = G.mono(sig, [0] * n)
            cur = P.add_move([(cur, Am)])
            # the term is D1 Am g ; its frame differs from D_ej g by the free left factor Am^-1 ... fine
    return P

def monomial_split(P):
    """every generator's coefficient vector as a polynomial in the parameters: list of (generator, monomial, vector).
    The span of all monomial parts (each its own generator) contains the C[B']-span of the program for EVERY choice of
    coefficients and every projective limit, so a target outside it at some f excludes the skeleton outright."""
    polys = []
    for i, T in enumerate(P.terms):
        if P.maps[i] is None:
            polys.append({(): np.ones(1, dtype=complex)}); continue
        out = {}
        for (pidx, src, si, ri) in P.maps[i]:
            for mono, vec in polys[src].items():
                m2 = tuple(sorted(mono + ((pidx,) if pidx is not None else ())))
                v = out.setdefault(m2, np.zeros(len(T), dtype=complex))
                np.add.at(v, ri, vec[si])
        polys.append(out)
    return polys

def split_test(P, fs=(2,), tol=1e-9):
    """residual of each output's target from the span of all monomial parts (diagonal blocks at f in fs)"""
    polys = monomial_split(P)
    res = {}
    for f in fs:
        d = 16 ** f
        rows = {}
        S = []
        for i, T in enumerate(P.terms):
            for mono, vec in polys[i].items():
                X = np.zeros((16, d * P.W), dtype=complex)
                for j, (l, g) in enumerate(T):
                    if abs(vec[j]) < 1e-14: continue
                    if g not in rows: rows[g] = np.array([DTest.tens(G.MATS[g][x], f) for x in range(16)])
                    X[:, d * l:d * (l + 1)] += vec[j] * rows[g]
                S.append(X)
        S = np.vstack(S)
        U, s, Vh = np.linalg.svd(S, full_matrices=False)
        Vr = Vh[: int((s > tol * s[0]).sum())]
        TK = np.array([DTest.tens(G.K[x], f) for x in range(16)])
        for l0 in range(P.W):
            Tm = np.zeros((16, d * P.W), dtype=complex); Tm[:, d * l0:d * (l0 + 1)] = TK
            Rm = Tm - (Tm @ Vr.conj().T) @ Vr
            res[(f, l0)] = round(float(np.sum(np.abs(Rm) ** 2) / np.sum(np.abs(Tm) ** 2)), 6)
    return res, sum(len(p) for p in polys)

class Closure:
    """FULL exact-structure test (bm.Game.closure, as mw2-absorb4/full_test.py): for output l0 the closure of the target
    cell over frames, with entries = term coefficients (functions of the parameters); residual = min_beta |A beta - e_K|^2.
    The structure is computed once with the supports of the skeleton (generic coefficients)."""
    def __init__(self, P, l0, cap=20000):
        self.P = P; self.l0 = l0
        tok = []; gd = []
        for i, T in enumerate(P.terms):
            dd = {}
            for j, tt in enumerate(T):
                dd[tt] = len(tok) + 1; tok.append((i, j))
            gd.append(dd)
        cl = G.closure(gd, (l0, G.KID), cap)
        self.ok = cl is not None
        if not self.ok: return
        rows, cols, ent = cl
        self.nr, self.nc = len(rows), len(cols)
        E = np.array([(r, c, tok[co - 1][0], tok[co - 1][1]) for (r, c, co) in ent], dtype=np.int64).reshape(-1, 4)
        sizes = [len(T) for T in P.terms]; off = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self.E = E; self.F = off[E[:, 2]] + E[:, 3]
        self.rows = rows; self.cols = cols
    def matrix(self, par):
        fl = np.concatenate(self.P.coefs(par))
        A = np.zeros((self.nr, self.nc), dtype=complex)
        np.add.at(A, (self.E[:, 0], self.E[:, 1]), fl[self.F])
        return A
    def resid_vec(self, par):
        A = self.matrix(par)
        t = np.zeros(self.nr, dtype=complex); t[0] = 1
        beta = np.linalg.lstsq(A, t, rcond=None)[0]
        return A @ beta - t
    def residual(self, par):
        return float(np.sum(np.abs(self.resid_vec(par)) ** 2))
