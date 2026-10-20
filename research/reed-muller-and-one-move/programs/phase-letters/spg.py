#!/usr/bin/env python3
"""spg.py (rm2-phase, 15 Oct 2026): Lemma N's block as a game on Sp(2m,2), with exact F_2 bookkeeping.

Vectors u in F_2^(2m) are ints: bits 0..m-1 = x, bits m..2m-1 = b (Pauli label X^x Z^b).  omega(u,v) = x.b' + b.x'.
A group element is the tuple of images of the 2m basis vectors.
  t_z   : (x|b) -> (x + (b.z) z | b)        the image of the paid move dir z (transvection along (z|0))
  x(S)  : (x|b) -> (x + S b | b)            so t_z = x(z z^T) and the target K = x(I) = prod_j t_(e_j)
  y(C)  : (x|b) -> (x | b + C x)            the images of the phase letters (C = w w^T for diag(i^[w.x]))
  lev(A): (x|b) -> (A x | A^-T b)           the images of the address maps
  w     : (x|b) -> (b|x)                    the long Weyl element
L0 = {(0|b)} (the frame of the identity: diagonal Paulis).  The free part is P = Stab(L0) = {lev(A) y(C)}.
"""
import itertools, random

def popc(v): return bin(v).count("1")

class Sp:
    def __init__(s, m):
        s.m = m; s.n = 2 * m; s.mask = (1 << m) - 1
    # --- vectors
    def X(s, u): return u & s.mask
    def B(s, u): return u >> s.m
    def vec(s, x, b): return x | (b << s.m)
    def om(s, u, v):
        return (popc(s.X(u) & s.B(v)) ^ popc(s.B(u) & s.X(v))) & 1
    # --- elements
    def from_fn(s, f): return tuple(f(1 << i) for i in range(s.n))
    def app(s, g, u):
        r = 0; i = 0
        while u:
            if u & 1: r ^= g[i]
            u >>= 1; i += 1
        return r
    def mul(s, g, h): return tuple(s.app(g, hi) for hi in h)      # (g h)(u) = g(h(u))
    def inv(s, g):
        # brute force over the group generated is too slow; invert as a linear map by Gaussian elimination
        n = s.n; rows = [(g[i], 1 << i) for i in range(n)]       # g(e_i) = g[i]; want g^-1(e_j)
        # solve: for each target e_j find combination of g[i] equal to e_j
        piv = {}
        basis = []
        for val, comb in rows:
            for (pv, pc) in basis:
                if val & (pv & -pv) and False: pass
            basis.append((val, comb))
        # standard elimination
        vals = [r[0] for r in rows]; combs = [r[1] for r in rows]
        for col in range(n):
            bit = 1 << col
            p = next((k for k in range(col, n) if vals[k] & bit), None)
            assert p is not None, "singular"
            vals[col], vals[p] = vals[p], vals[col]; combs[col], combs[p] = combs[p], combs[col]
            for k in range(n):
                if k != col and vals[k] & bit:
                    vals[k] ^= vals[col]; combs[k] ^= combs[col]
        # now vals[col] = e_col and combs[col] = sum of e_i with g(sum) = e_col
        return tuple(combs[j] for j in range(n))
    def ident(s): return tuple(1 << i for i in range(s.n))
    def is_symplectic(s, g):
        return all(s.om(g[i], g[j]) == s.om(1 << i, 1 << j) for i in range(s.n) for j in range(s.n))
    def matvec(s, M, x):           # M: list of m row-ints (row i = bits of row i); returns M x
        return sum(((popc(M[i] & x) & 1) << i) for i in range(s.m))
    def t(s, z):
        return s.from_fn(lambda u: u ^ (s.vec(z, 0) if popc(s.B(u) & z) & 1 else 0))
    def xS(s, S): return s.from_fn(lambda u: s.vec(s.X(u) ^ s.matvec(S, s.B(u)), s.B(u)))
    def yC(s, C): return s.from_fn(lambda u: s.vec(s.X(u), s.B(u) ^ s.matvec(C, s.X(u))))
    def lev(s, A):
        Ainv = inv_f2(A, s.m); AinvT = transpose(Ainv, s.m)
        return s.from_fn(lambda u: s.vec(s.matvec(A, s.X(u)), s.matvec(AinvT, s.B(u))))
    def w(s): return s.from_fn(lambda u: s.vec(s.B(u), s.X(u)))
    def I(s): return [1 << i for i in range(s.m)]
    def K(s): return s.xS(s.I())
    # --- frames
    def span(s, vs):
        S = {0}
        for v in vs:
            if v not in S: S |= {a ^ v for a in S}
        return frozenset(S)
    def L0(s): return s.span([s.vec(0, 1 << i) for i in range(s.m)])
    def LX(s): return s.span([s.vec(1 << i, 0) for i in range(s.m)])
    def act(s, g, U): return frozenset(s.app(g, u) for u in U)
    def dist(s, U, V):
        k = len(U & V).bit_length() - 1
        return s.m - k

def transpose(M, m):
    return [sum((((M[r] >> c) & 1) << r) for r in range(m)) for c in range(m)]

def inv_f2(M, m):
    rows = [M[i] | (1 << (m + i)) for i in range(m)]
    for c in range(m):
        p = next(k for k in range(c, m) if (rows[k] >> c) & 1)
        rows[c], rows[p] = rows[p], rows[c]
        for k in range(m):
            if k != c and (rows[k] >> c) & 1: rows[k] ^= rows[c]
    return [r >> m for r in rows]

def rand_gl(m, rng):
    while True:
        M = [rng.randrange(1 << m) for _ in range(m)]
        try:
            inv_f2(M, m); return M
        except StopIteration:
            pass

def rand_sym(m, rng):
    S = [0] * m
    for i in range(m):
        for j in range(i, m):
            if rng.random() < 0.5:
                S[i] |= 1 << j; S[j] |= 1 << i
    return S

def lagrangians(sp):
    """all Lagrangians, by BFS from L0 under the generators t_(e_i), lev(elementary), y(e_i e_i^T)."""
    gens = [sp.t(1 << i) for i in range(sp.m)]
    gens += [sp.yC([(1 << i) if r == i else 0 for r in range(sp.m)]) for i in range(sp.m)]
    m = sp.m
    for i in range(m):
        for j in range(m):
            if i != j:
                E = [(1 << r) | ((1 << j) if r == i else 0) for r in range(m)]
                gens.append(sp.lev(E))
    start = sp.L0(); seen = {start: 0}; order = [start]
    k = 0
    while k < len(order):
        U = order[k]; k += 1
        for g in gens:
            V = sp.act(g, U)
            if V not in seen: seen[V] = len(order); order.append(V)
    return order, seen
