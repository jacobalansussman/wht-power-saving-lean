#!/usr/bin/env python3
"""frames.py (p6-prove, 15 Oct 2026): the frame graph Gamma_m (Lagrangians of F_2^{2m}) and Sp(2m,2) acting on it.

Vectors v = (a|b) of F_2^{2m} are ints: a = low m bits, b = high m bits.  omega(v,w) = a_v.b_w + a_w.b_v (mod 2).
L0 = {(0|b)} (raw frame), LT = {(z|z)} (frame of C^{(x)m}); d(A,B) = m - dim(A cap B).
An Sp element is a tuple of the images of the 2m basis vectors.  Exact (pure integer) code.
"""
import itertools, random
from collections import deque

def pc(x): return bin(x).count("1") & 1

class Frames:
    def __init__(s, m):
        s.m = m; s.n = 2 * m; s.mask = (1 << m) - 1
        L0 = s.span([1 << (m + i) for i in range(m)])
        LT = s.span([(1 << i) | (1 << (m + i)) for i in range(m)])
        # all Lagrangians: orbit of L0 under transvections T_v (v ranges over all nonzero vectors)
        s.vecs = list(range(1, 1 << s.n))
        s.pts = [L0]; s.idx = {L0: 0}
        q = deque([L0])
        while q:
            A = q.popleft()
            for v in s.vecs:
                B = frozenset(s.transvect(v, x) for x in A)
                if B not in s.idx: s.idx[B] = len(s.pts); s.pts.append(B); q.append(B)
        s.N = len(s.pts)
        s.L0 = 0; s.LT = s.idx[LT]
        # adjacency: dim(A cap B) = m - 1  <=> |A cap B| = 2^(m-1)
        s.nb = [[] for _ in range(s.N)]
        half = 1 << (m - 1)
        if s.N <= 200:
            for i in range(s.N):
                for j in range(i + 1, s.N):
                    if len(s.pts[i] & s.pts[j]) == half: s.nb[i].append(j); s.nb[j].append(i)
        else:
            # neighbours of A are the images of A under transvections T_v with v not in A (they meet A in a hyperplane)
            for i, A in enumerate(s.pts):
                nbs = set()
                for v in s.vecs:
                    if v in A: continue
                    B = frozenset(s.transvect(v, x) for x in A)
                    j = s.idx[B]
                    if j != i: nbs.add(j)
                s.nb[i] = sorted(nbs)
        # all-pairs distances only for small graphs; otherwise distances from L0 (and LT) only; d() is exact anyway
        s.dist = [s.bfs(i) for i in range(s.N)] if s.N <= 200 else {0: s.bfs(0), s.LT: s.bfs(s.LT)}
    def omega(s, v, w):
        m = s.mask
        return pc((v & m) & (w >> s.m)) ^ pc((w & m) & (v >> s.m))
    def transvect(s, v, x): return x ^ v if s.omega(x, v) else x
    def span(s, gens):
        S = {0}
        for g in gens: S |= {x ^ g for x in S}
        return frozenset(S)
    def bfs(s, i):
        d = [-1] * s.N; d[i] = 0; q = deque([i])
        while q:
            a = q.popleft()
            for b in s.nb[a]:
                if d[b] < 0: d[b] = d[a] + 1; q.append(b)
        return d
    def lv(s, i): return s.dist[0][i]
    def d(s, i, j): return s.m - (len(s.pts[i] & s.pts[j]).bit_length() - 1)
    # ---- Sp elements as images of basis vectors
    def apply_vec(s, g, x):
        y = 0; k = 0
        while x:
            if x & 1: y ^= g[k]
            x >>= 1; k += 1
        return y
    def apply_pt(s, g, i): return s.idx[frozenset(s.apply_vec(g, x) for x in s.pts[i])]
    def mul(s, g, h): return tuple(s.apply_vec(g, h[k]) for k in range(s.n))      # g after h
    def ident(s): return tuple(1 << k for k in range(s.n))
    def tv(s, v): return tuple(s.transvect(v, 1 << k) for k in range(s.n))
    def is_symp(s, g):
        return all(s.omega(g[i], g[j]) == s.omega(1 << i, 1 << j) for i in range(s.n) for j in range(s.n))
    def perm(s, g): return tuple(s.apply_pt(g, i) for i in range(s.N))
    def all_elements(s):
        """whole group as permutations of the points (only for m <= 2)"""
        gens = [s.perm(s.tv(v)) for v in s.vecs]
        e = tuple(range(s.N)); seen = {e}; q = deque([e])
        while q:
            p = q.popleft()
            for g in gens:
                r = tuple(g[p[i]] for i in range(s.N))
                if r not in seen: seen.add(r); q.append(r)
        return sorted(seen)
    def rand_elem(s, rng, steps=40):
        g = s.ident()
        for _ in range(steps): g = s.mul(s.tv(rng.choice(s.vecs)), g)
        return g

if __name__ == "__main__":
    for m in (1, 2, 3):
        F = Frames(m)
        deg = {len(x) for x in F.nb}
        print("m=%d points=%d degree=%s d(L0,LT)=%d diameter=%d" % (m, F.N, deg, F.dist[0][F.LT], max(F.dist[0])))
        if m <= 2:
            G = F.all_elements(); print("   |Sp| as permutations of points:", len(G))
        rng = random.Random(1)
        g = F.rand_elem(rng); assert F.is_symp(g)
        # isometry check on random element
        p = F.perm(g)
        ok = all(F.dist[i][j] == F.dist[p[i]][p[j]] for i in range(F.N) for j in range(0, F.N, 7))
        print("   random element symplectic, isometry:", ok)
