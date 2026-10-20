#!/usr/bin/env python3
"""veng.py (v-p6-prove, 15 Oct 2026): an independent exact engine for the engine game at block m, TRUE model.

Written from scratch (no import of ledger.py or p6-prove code).  Group elements are exact integer matrices
g * 2^m (re, im int64), keyed WITHOUT normalising scalars: a scalar matrix produced by letters (e.g. XZXZ = -I)
is a group element of its own, distinct from the coefficient -1 (the every-f model C[G], sl-alg Fact 0).

Letters (sl-comb PROOF.md section 0):
  paid  D<z>      dir z = (1+i)/2 * 1 + (1-i)/2 * X^z, z != 0
  free  X<z>      shift x -> x ^ z
        A<c0,..>  invertible address map, column images
        P<z>_<c>  diag(i^((z.x + c) mod 2))
        S<k>      scalar matrix i^k * I  (a group element)
Steps: ("move", r, z), ("free", r, name), ("mix", r, [(s, mu, name|None), ...]) (old values; ledger semantics).
Contents: dict {(input l, gkey): coefficient (complex Fraction pair)}.
Frames: position lam(g) = gbar^-1 L0, L0 = span{(0|e_i)}; computed from my own Pauli-conjugation code.
"""
import numpy as np
from fractions import Fraction as Fr

class Q:
    """Gaussian rational"""
    __slots__ = ("a", "b")
    def __init__(s, a, b=0): s.a = Fr(a); s.b = Fr(b)
    def __add__(s, o): return Q(s.a + o.a, s.b + o.b)
    def __mul__(s, o): return Q(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)
    def z(s): return s.a == 0 and s.b == 0
    def abs2(s): return s.a * s.a + s.b * s.b
    def c(s): return complex(float(s.a), float(s.b))
    def __repr__(s): return "(%s%+si)" % (s.a, s.b)

def toQ(x):
    if isinstance(x, Q): return x
    if isinstance(x, complex): return Q(Fr(x.real), Fr(x.imag))
    return Q(x)

class Eng:
    def __init__(s, m):
        s.m = m; s.N = 1 << m; s.S = 1 << m
        s.mats = {}          # key -> (re, im)
        s.cache = {}
        s.symcache = {}
        N, S = s.N, s.S
        I = (np.eye(N, dtype=np.int64) * S, np.zeros((N, N), dtype=np.int64))
        s.ID = s.key(I)
        s.letters = {}
        for z in range(1, N):
            s.letters["X%d" % z] = s.key(s.perm(lambda x, z=z: x ^ z))
            re = np.zeros((N, N), dtype=np.int64); im = np.zeros((N, N), dtype=np.int64)
            h = S // 2
            for x in range(N):
                re[x, x] += h; im[x, x] += h
                re[x ^ z, x] += h; im[x ^ z, x] -= h
            s.letters["D%d" % z] = s.key((re, im))
        for z in range(N):
            for c in (0, 1):
                re = np.zeros((N, N), dtype=np.int64); im = np.zeros((N, N), dtype=np.int64)
                for x in range(N):
                    if (bin(z & x).count("1") + c) % 2: im[x, x] = S
                    else: re[x, x] = S
                s.letters["P%d_%d" % (z, c)] = s.key((re, im))
        for k in range(4):
            u = [(1, 0), (0, 1), (-1, 0), (0, -1)][k]
            s.letters["S%d" % k] = s.key((np.eye(N, dtype=np.int64) * S * u[0], np.eye(N, dtype=np.int64) * S * u[1]))
        g = s.ID
        for j in range(m): g = s.lmul(s.letters["D%d" % (1 << j)], g)
        s.K = g
    def perm(s, f):
        re = np.zeros((s.N, s.N), dtype=np.int64)
        for x in range(s.N): re[f(x), x] = s.S
        return (re, np.zeros((s.N, s.N), dtype=np.int64))
    def addr(s, cols):
        name = "A" + ",".join(map(str, cols))
        if name not in s.letters:
            def f(x):
                y = 0
                for k in range(s.m):
                    if (x >> k) & 1: y ^= cols[k]
                return y
            assert len({f(x) for x in range(s.N)}) == s.N
            s.letters[name] = s.key(s.perm(f))
        return name
    def key(s, M):
        k = M[0].tobytes() + M[1].tobytes()
        if k not in s.mats: s.mats[k] = (M[0].copy(), M[1].copy())
        return k
    def mulk(s, a, b):
        ck = (a, b)
        if ck not in s.cache:
            A, B = s.mats[a], s.mats[b]
            re = A[0] @ B[0] - A[1] @ B[1]; im = A[0] @ B[1] + A[1] @ B[0]
            assert (re % s.S == 0).all() and (im % s.S == 0).all()
            s.cache[ck] = s.key((re // s.S, im // s.S))
        return s.cache[ck]
    def lmul(s, a, b): return s.mulk(a, b)
    def dag(s, a):
        A = s.mats[a]; return s.key((A[0].T.copy(), -A[1].T.copy()))
    def cplx(s, a):
        A = s.mats[a]; return (A[0] + 1j * A[1]) / s.S
    # ---------------- symplectic image, by conjugating Paulis (integers only)
    def pauli(s, a, b):
        re = np.zeros((s.N, s.N), dtype=np.int64)
        for x in range(s.N): re[x ^ a, x] = -1 if bin(b & x).count("1") % 2 else 1
        return (re * s.S, np.zeros((s.N, s.N), dtype=np.int64))
    def symp(s, g):
        """columns = images of basis labels; label u = a | (b << m), Pauli X^a Z^b; returns list of 2m ints"""
        if g in s.symcache: return s.symcache[g]
        m, N = s.m, s.N
        gd = s.dag(g)
        cols = []
        for j in range(2 * m):
            a, b = (1 << j, 0) if j < m else (0, 1 << (j - m))
            Pk = s.key(s.pauli(a, b))
            Qk = s.mulk(s.mulk(g, Pk), gd)
            re, im = s.mats[Qk]
            nzc = [(x, y) for x in range(N) for y in range(N) if re[x, y] or im[x, y]]
            assert len(nzc) == N, "not monomial: not Clifford"
            a2 = [x for (x, y) in nzc if y == 0][0]
            ph = (re[a2, 0], im[a2, 0])
            b2 = 0
            for k in range(m):
                xk = 1 << k
                v = (re[xk ^ a2, xk], im[xk ^ a2, xk])
                if v == (-ph[0], -ph[1]): b2 |= xk
                else: assert v == ph
            cols.append(a2 | (b2 << m))
        s.symcache[g] = cols
        return cols
    def apply_symp(s, cols, v):
        out = 0
        for j in range(2 * s.m):
            if (v >> j) & 1: out ^= cols[j]
        return out
    def span(s, vecs):
        """canonical fully reduced echelon basis (sorted tuple) of the F2 span of 2m-bit ints"""
        basis = {}
        for v in vecs:
            while v:
                p = v.bit_length() - 1
                if p in basis: v ^= basis[p]
                else: basis[p] = v; break
        for p in sorted(basis, reverse=True):
            for q in basis:
                if q != p and (basis[q] >> p) & 1: basis[q] ^= basis[p]
        return tuple(sorted(basis.values()))
    def L0(s): return s.span([1 << (s.m + i) for i in range(s.m)])
    def pos(s, g):
        """lam(g) = gbar^-1 L0"""
        cols = s.symp(s.dag(g))
        return s.span([s.apply_symp(cols, 1 << (s.m + i)) for i in range(s.m)])
    def rank(s, vecs): return len(s.span(list(vecs)))
    def d(s, A, B): return s.rank(list(A) + list(B)) - s.m
    def level(s, g): return s.d(s.L0(), s.pos(g))

def step(E, C, st):
    C = [dict(c) for c in C]
    kind = st[0]
    def act(gname, c):
        out = {}
        for (l, g), co in c.items():
            k = (l, E.mulk(gname, g)); out[k] = out.get(k, Q(0)) + co
        return {k: v for k, v in out.items() if not v.z()}
    if kind == "move":
        C[st[1]] = act(E.letters["D%d" % st[2]], C[st[1]])
    elif kind == "free":
        C[st[1]] = act(E.letters[st[2]], C[st[1]])
    elif kind == "mix":
        new = {}
        for (sr, mu, name) in st[2]:
            mu = toQ(mu)
            src = C[sr] if name is None else act(E.letters[name], C[sr])
            for k, co in src.items(): new[k] = new.get(k, Q(0)) + mu * co
        C[st[1]] = {k: v for k, v in new.items() if not v.z()}
    return C

def run(E, W, R, word, deliver=None):
    """returns (is_unit_true_model, M, history of contents after each step)"""
    C = [dict() for _ in range(R)]
    for r in range(W): C[r] = {(r, E.ID): Q(1)}
    hist = [C]
    for st in word:
        C = step(E, C, st); hist.append(C)
    M = sum(1 for st in word if st[0] == "move")
    deliver = deliver or list(range(W))
    scal = {E.letters["S%d" % k] for k in range(4)}
    targets = {E.mulk(sk, E.K) for sk in scal}
    ok = True
    for l, r in enumerate(deliver):
        c = C[r]
        if len(c) != 1: ok = False; continue
        (ll, g), co = next(iter(c.items()))
        if ll != l or g not in targets or co.abs2() != 1: ok = False
    return ok, M, hist

def raw_check(E, W, R, word, f, rng, deliver=None):
    """raw complex arrays at f digit columns: letters act as g^(x)f; returns max relative error vs u*K^(x)f x_l, |u| = 1"""
    n = E.N
    shape = (n,) * f
    def apply(gk, v):
        g = E.cplx(gk)
        for ax in range(f):
            v = np.moveaxis(np.tensordot(g, v, axes=([1], [ax])), 0, ax)
        return v
    xs = [rng.standard_normal(shape) + 1j * rng.standard_normal(shape) for _ in range(W)]
    regs = [xs[r].copy() if r < W else np.zeros(shape, dtype=complex) for r in range(R)]
    for st in word:
        if st[0] == "move": regs[st[1]] = apply(E.letters["D%d" % st[2]], regs[st[1]])
        elif st[0] == "free": regs[st[1]] = apply(E.letters[st[2]], regs[st[1]])
        else:
            new = np.zeros(shape, dtype=complex)
            for (sr, mu, name) in st[2]:
                src = regs[sr] if name is None else apply(E.letters[name], regs[sr])
                new = new + toQ(mu).c() * src
            regs[st[1]] = new
    deliver = deliver or list(range(W))
    worst = 0.0
    for l, r in enumerate(deliver):
        t = apply(E.K, xs[l])
        u = np.vdot(t, regs[r]) / np.vdot(t, t)
        err = np.linalg.norm(regs[r] - u * t) / np.linalg.norm(t)
        worst = max(worst, err, abs(abs(u) - 1))
    return worst
