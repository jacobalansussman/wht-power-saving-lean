#!/usr/bin/env python3
"""sp.py: Sp(2m,2) on Lagrangians (frames).  Vectors of F_2^{2m} are ints: bits 0..m-1 = x, bits m..2m-1 = b.
omega((x,b),(x',b')) = x.b' + b.x'.  L0 = {(0,b)} (the frame the free group B stabilises), distance of a frame L from
L0 = m - dim(L cap L0).  Moves: transvection t = T_v, v = (e_1, 0).  Target Kbar: (x,b) -> (x+b, b).
Checks: number of frames = prod(2^j+1); the free group P = Stab(L0) has exactly m+1 orbits on frames (one per distance)
=> Hecke algebra e C[Sp] e has dimension m+1 (Gelfand pair, multiplicity-free frame module); the dual polar graph's
intersection array b_i = 2^{i+1}(2^{m-i}-1), c_i = 2^i - 1; spherical values omega_j(d); abelianisation of P."""
import itertools, sys
import numpy as np

def pc(v): return bin(v).count("1") & 1
class Sp:
    def __init__(s, m):
        s.m = m; s.n = 2 * m; s.mx = (1 << m) - 1
    def om(s, u, v):
        m = s.m
        return pc((u & s.mx) & (v >> m)) ^ pc((u >> m) & (v & s.mx))
    def span(s, gens):
        S = {0}
        for g in gens: S |= {x ^ g for x in S}
        return frozenset(S)
    def apply(s, M, v):  # M: list of images of basis vectors (ints)
        out = 0
        for j in range(s.n):
            if (v >> j) & 1: out ^= M[j]
        return out
    def mat_mul(s, A, B):  # A after B
        return [s.apply(A, B[j]) for j in range(s.n)]
    def transvection(s, v):
        return [(1 << j) ^ (v if s.om(1 << j, v) else 0) for j in range(s.n)]
    def frames(s):
        L0 = s.span([1 << (s.m + j) for j in range(s.m)])
        gens = [s.transvection(v) for v in range(1, 1 << s.n)]
        seen = {L0}; todo = [L0]
        while todo:
            L = todo.pop()
            for g in gens:
                L2 = frozenset(s.apply(g, x) for x in L)
                if L2 not in seen: seen.add(L2); todo.append(L2)
        return L0, sorted(seen, key=lambda L: sorted(L))
    def dist(s, L, L0):
        return s.m - (len(L & L0).bit_length() - 1)
    def P_gens(s):
        """generators of P = Stab(L0): GL(m,2) acting as (x,b) -> (Ax, A^{-T} b) (elementary A) and shears
        (x,b) -> (x, b + S x) with S symmetric (elementary S)."""
        m = s.m; gens = []
        for i in range(m):
            for j in range(m):
                if i == j: continue
                # A = I + E_ij : x_i += x_j ; A^{-T} = I + E_ji : b_j += b_i
                M = []
                for k in range(s.n):
                    v = 1 << k
                    x, b = v & s.mx, v >> m
                    x2 = x ^ ((((x >> j) & 1) << i))
                    b2 = b ^ ((((b >> i) & 1) << j))
                    M.append(x2 | (b2 << m))
                gens.append(M)
        for i in range(m):
            for j in range(i, m):
                M = []
                for k in range(s.n):
                    v = 1 << k
                    x, b = v & s.mx, v >> m
                    # b += S x with S = E_ij + E_ji (or E_ii)
                    add = 0
                    if (x >> j) & 1: add ^= 1 << i
                    if i != j and (x >> i) & 1: add ^= 1 << j
                    M.append(x | ((b ^ add) << m))
                gens.append(M)
        return gens

def main(m):
    S = Sp(m)
    L0, F = S.frames()
    idx = {L: i for i, L in enumerate(F)}
    n = len(F)
    expect = 1
    for j in range(1, m + 1): expect *= (2 ** j + 1)
    # check P generators stabilise L0 and preserve omega
    Pg = S.P_gens()
    for g in Pg:
        assert frozenset(S.apply(g, x) for x in L0) == L0
        for u in range(1 << S.n):
            for v in range(1 << S.n):
                if S.om(S.apply(g, u), S.apply(g, v)) != S.om(u, v): raise SystemExit("not symplectic")
    # P-orbits on frames
    perm = [[idx[frozenset(S.apply(g, x) for x in L)] for L in F] for g in Pg]
    orb = [-1] * n; no = 0
    for i0 in range(n):
        if orb[i0] >= 0: continue
        orb[i0] = no; st = [i0]
        while st:
            i = st.pop()
            for p in perm:
                if orb[p[i]] < 0: orb[p[i]] = no; st.append(p[i])
        no += 1
    dists = [S.dist(L, L0) for L in F]
    orbit_is_distance = all((orb[i] == orb[j]) == (dists[i] == dists[j]) for i in range(n) for j in range(0, n, max(1, n // 50)))
    print("m=%d: frames %d (expected %d); P-orbits on frames: %d (m+1 = %d); orbits = distance classes: %s; sizes %s"
          % (m, n, expect, no, m + 1, orbit_is_distance, [dists.count(d) for d in range(m + 1)]))
    # intersection array from the graph at the base frame
    adj = [[S.dist(F[i], F[j]) == 1 for j in range(n)] for i in range(n)] if n <= 2300 else None
    if adj is not None:
        i0 = idx[L0]
        b = []; c = []
        for d in range(m + 1):
            Ld = [j for j in range(n) if dists[j] == d]
            j = Ld[0]
            nb = [k for k in range(n) if adj[j][k]]
            c.append(sum(1 for k in nb if dists[k] == d - 1)); b.append(sum(1 for k in nb if dists[k] == d + 1))
        print("   intersection array b =", b, " c =", c, " formula b_i = 2^{i+1}(2^{m-i}-1):",
              [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)], " c_i = 2^i - 1:", [2 ** i - 1 for i in range(m + 1)])
    # transvection t and Kbar distances
    t = S.transvection(1)  # v = (e_1, 0)
    Kb = [((1 << j) if j < m else ((1 << j) ^ (1 << (j - m)))) for j in range(S.n)]
    for name, g in (("t", t), ("Kbar", Kb)):
        L = frozenset(S.apply(g, x) for x in L0)
        print("   %s maps L0 to a frame at distance %d" % (name, S.dist(L, L0)))
    # slice comparison: characters of the slice group trivial on shifts = binary quadratic forms = 2^{m(m+1)/2}
    print("   slice: the trivial shift character sits in 2^(m(m+1)/2) = %d characters; engine: the trivial P-component sits in %d irreps"
          % (2 ** (m * (m + 1) // 2), no))
    return S, F, L0, dists, perm

def spherical(m):
    """spherical values omega_j(d) of the dual polar graph of Sp(2m,2), from the intersection array."""
    b = [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)]
    c = [2 ** i - 1 for i in range(m + 1)]
    k = b[0]
    a = [k - b[i] - c[i] for i in range(m + 1)]
    L = np.zeros((m + 1, m + 1))
    for i in range(m + 1):
        L[i, i] = a[i]
        if i > 0: L[i, i - 1] = c[i]
        if i < m: L[i, i + 1] = b[i]
    th = np.sort(np.linalg.eigvals(L).real)[::-1]
    om = []
    for t_ in th:
        w = [1.0, t_ / k]
        for i in range(1, m):
            w.append((t_ * w[i] - c[i] * w[i - 1] - a[i] * w[i]) / b[i])
        om.append(w)
    return th, np.array(om)

if __name__ == "__main__":
    for m in (2, 3, 4):
        main(m)
    for m in (2, 3, 4, 5, 6):
        th, om = spherical(m)
        print("m=%d eigenvalues %s ; omega_j(1) = %s ; omega_j(m) = %s" % (m, np.round(th, 6), np.round(om[:, 1], 6), np.round(om[:, m], 6)))
