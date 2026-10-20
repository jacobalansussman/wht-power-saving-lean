#!/usr/bin/env python3
"""vlem43.py (v-p6-refute): check Lemma 4.3 of p6-refute with actual symplectic matrices over F_2.
Space F_2^(2m), basis e_1..e_m, f_1..f_m, omega(e_i, f_j) = delta_ij; L0 = span f (raw frame), L_T = span(e_i + f_i).
Cube B_m: hyperoctahedral permutations of 2m points, point 2i <-> f_i, 2i+1 <-> e_i.
Hamming H(m,3): S_3^m on points 3i <-> f_i, 3i+1 <-> e_i, 3i+2 <-> e_i + f_i (the three lines of plane i).
Checks: every group element is a symplectic matrix (homomorphism); the sub-game's free group = sub-group cap Stab(L0);
position P g -> frame g^-1 L0 is well defined and injective; the frame-graph distance m - dim(A cap B) equals the
coset-graph distance for every pair of positions; every paid letter h has d(L0, h L0) = 1 (distance-1 double coset of
Stab(L0) in Sp(2m,2)); the twist action P g -> P g t equals frame -> t^-1 frame; the sub-game's level-m frames are
at frame distance m from L0 (so they lie in the double coset of the transform)."""
import sys
import vgroups as VG

def rank2(vs):
    vs = list(vs); r = 0; basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v: basis.append(v); r += 1
    return r

def span(vs):
    out = {0}
    for v in vs:
        out |= {x ^ v for x in out}
    return frozenset(out)

def check(kind, m):
    n2 = 2 * m
    E = lambda i: 1 << i; F = lambda i: 1 << (m + i)            # e_i bit i, f_i bit m+i
    def omega(u, v):
        return (bin((u & ((1 << m) - 1)) & (v >> m)).count("1") + bin((u >> m) & (v & ((1 << m) - 1))).count("1")) & 1
    if kind == "cube":
        g = VG.hyperoct(m)
        vec = {}
        for i in range(m): vec[2 * i] = F(i); vec[2 * i + 1] = E(i)
    else:
        g = VG.hamming(m)
        vec = {}
        for i in range(m): vec[3 * i] = F(i); vec[3 * i + 1] = E(i); vec[3 * i + 2] = E(i) | F(i)
    def mat(p):
        # linear map sending e_i -> vec[p[pt(e_i)]], f_i -> vec[p[pt(f_i)]]
        inv = {v: k for k, v in vec.items()}
        img = {}
        for i in range(m):
            img[E(i)] = vec[p[inv[E(i)]]]; img[F(i)] = vec[p[inv[F(i)]]]
        def ap(v):
            r = 0
            for b in range(n2):
                if (v >> b) & 1: r ^= img[1 << b]
            return r
        return ap
    maps = [mat(p) for p in g.elems]
    bad = 0
    basis = [E(i) for i in range(m)] + [F(i) for i in range(m)]
    # homomorphism on the point action (for the Hamming case: the map must send every line point consistently)
    for k, p in enumerate(g.elems):
        ap = maps[k]
        for pt, v in vec.items():
            if ap(v) != vec[p[pt]]: bad += 1; break
        for u in basis:
            for v in basis:
                if omega(ap(u), ap(v)) != omega(u, v): bad += 1
    L0 = span([F(i) for i in range(m)])
    LT = span([E(i) | F(i) for i in range(m)])
    def frame(e):   # g^-1 L0
        ap = maps[g.I[e]]
        return frozenset(ap(v) for v in L0)
    stab = [e for e in range(g.n) if frame(e) == L0]
    if sorted(stab) != sorted(g.P): bad += 1; print("free group != subgroup cap Stab(L0)")
    fr = {}
    for e in range(g.n):
        c = g.pos[e]
        F_ = frame(e)
        if c in fr and fr[c] != F_: bad += 1
        fr[c] = F_
    if len(set(fr.values())) != g.N: bad += 1; print("not injective")
    d = lambda A, B: m - rank2(A & B)
    # coset-graph distances
    import collections
    for a in range(g.N):
        dist = {a: 0}; q = collections.deque([a])
        while q:
            x = q.popleft()
            for y in g.nb[x]:
                if y not in dist: dist[y] = dist[x] + 1; q.append(y)
        for b in range(g.N):
            if dist[b] != d(fr[a], fr[b]): bad += 1
    for h in g.H:
        if d(L0, frozenset(maps[h](v) for v in L0)) != 1: bad += 1
    # twist action
    for t in range(g.n):
        ti = maps[g.I[t]]
        for c in range(g.N):
            if fr[g.ract[t][c]] != frozenset(ti(v) for v in fr[c]): bad += 1
    tops = [c for c in range(g.N) if g.lv[c] == g.m]
    for c in tops:
        if d(L0, fr[c]) != m: bad += 1
    print("%s m=%d: |G| %d, positions %d, paid letters %d, level-m frames %d (L_T among them: %s), failures %d" % (
        kind, m, g.n, g.N, len(g.H), len(tops), any(fr[c] == LT for c in tops), bad))

if __name__ == "__main__":
    for kind, m in (("cube", 2), ("cube", 3), ("cube", 4), ("hamming", 2), ("hamming", 3), ("hamming", 4)):
        check(kind, m)
