#!/usr/bin/env python3
"""sp.py (w8-rmexact): Lemma N's Pauli-trivial block C[Sp(2m,2)] - deletion and orbit-union Theorem R on permutation
modules.  V = F_2^{2m}, u = (x|b) (x bits 0..m-1, b bits m..2m-1), omega(u,v) = x.b' + b.x'.  Move: the transvection
t_z (x|b) -> (x + (b.z) z | b) (vector (z|0)); target K (x|b) -> (x + b | b) = prod_j t_{e_j}; free part: the Siegel
parabolic P_{L0}, L0 = {(0|b)}... (Lemma N: stabiliser of the Z-type Lagrangian; here L0 = {(x|0)}? we take the
stabiliser of the Lagrangian Z = {(0|b)} AND of X = {(x|0)} as two options and report both).
For each permutation module (G-set: nonzero vectors, isotropic k-spaces, ...), decompose into irreducible constituents
(eigenspaces of a generic element of the orbital algebra), and print per constituent: dim, r_t = dim of the -1
eigenspace of t, r_K, ratio r_K / r_t (deletion: M r_t >= W r_K).  Theorem R: Y_0 = indicators of unions of
P-orbits on the G-set, f(g) = |gU \\ U|."""
import sys, itertools, numpy as np
m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
which = sys.argv[2].split(',') if len(sys.argv) > 2 else ['pts', 'lag', 'iso2']
N = 2 * m; nz = list(range(1, 1 << N))
def om(u, v): return bin(((u & ((1 << m) - 1)) & (v >> m)) ^ ((u >> m) & (v & ((1 << m) - 1)))).count('1') & 1
def trans(v): return lambda u: u ^ v if om(u, v) else u
t = trans(1)                                    # t_{e1}: vector (e1|0)
def K(u): return u ^ (u >> m)                   # (x|b) -> (x+b|b)
gens = [trans(v) for v in nz]                   # all transvections generate Sp(2m,2)
def span(vs):
    s = {0}
    for v in vs: s |= {a ^ v for a in s}
    return frozenset(s)
def iso_spaces(k):
    out = set()
    def rec(basis, start):
        if len(basis) == k: out.add(span(basis)); return
        S = span(basis)
        for v in range(start, 1 << N):
            if v in S: continue
            if all(om(v, w) == 0 for w in basis): rec(basis + [v], v + 1)
    rec([], 1)
    return sorted(out, key=lambda s: sorted(s))
def module(name):
    if name == 'pts': X = [frozenset([0, v]) for v in nz]
    elif name == 'lag': X = iso_spaces(m)
    elif name.startswith('iso'): X = iso_spaces(int(name[3:]))
    pos = {s: i for i, s in enumerate(X)}
    act = lambda g, s: frozenset(g(v) for v in s)
    P = [[pos[act(g, s)] for s in X] for g in gens]
    pt = [pos[act(t, s)] for s in X]; pK = [pos[act(K, s)] for s in X]
    return X, P, pt, pK, act, pos
for name in which:
    X, P, pt, pK, act, pos = module(name); n = len(X)
    # orbitals by BFS on pairs
    lab = -np.ones((n, n), dtype=np.int64); no = 0
    for a in range(n):
        for b in range(n):
            if lab[a, b] >= 0: continue
            st = [(a, b)]; lab[a, b] = no
            while st:
                i, j = st.pop()
                for p in P:
                    i2, j2 = p[i], p[j]
                    if lab[i2, j2] < 0: lab[i2, j2] = no; st.append((i2, j2))
            no += 1
    rng = np.random.default_rng(3)
    r = rng.normal(size=no) + 1j * rng.normal(size=no)
    T = r[lab]; T = T + T.conj().T
    w, V = np.linalg.eigh(T)
    groups = []; cur = [0]
    for i in range(1, n):
        if abs(w[i] - w[i - 1]) < 1e-7 * max(1, abs(w[i])): cur.append(i)
        else: groups.append(cur); cur = [i]
    groups.append(cur)
    Pt = np.zeros((n, n)); Pt[pt, range(n)] = 1
    PK = np.zeros((n, n)); PK[pK, range(n)] = 1
    print("module %s: %d points, rank %d, constituents %d (multiplicity-free %s)" % (name, n, no, len(groups), no == len(groups)))
    best = 0
    for g in groups:
        B = V[:, g]; d = len(g)
        trt = np.trace(B.conj().T @ Pt @ B).real; trK = np.trace(B.conj().T @ PK @ B).real
        rt = (d - trt) / 2; rK = (d - trK) / 2
        ratio = rK / rt if rt > 1e-9 else (float('inf') if rK > 1e-9 else 0)
        best = max(best, ratio)
        print("   dim %4d  chi(t) %7.2f  chi(K) %7.2f  r_t %6.1f  r_K %6.1f  ratio %.4f" % (d, trt, trK, rt, rK, ratio))
    print("   best deletion ratio in this module: %.4f" % best)
    # Theorem R with orbit unions of the two parabolics
    for pname, L in (('stab Z={(0|b)}', [v << m for v in range(1, 1 << m)]), ('stab X={(x|0)}', list(range(1, 1 << m)))):
        Lset = span(L)
        Pg = [p for p, g in zip(P, gens) if all(g(v) in Lset for v in L)]
        # transvections preserving L generate only part of the parabolic; add GL action on L too
        Lev = []
        for perm in itertools.permutations(range(m)):
            def g(u, perm=perm):
                x = u & ((1 << m) - 1); b = u >> m
                px = sum(((x >> i) & 1) << perm[i] for i in range(m)); pb = sum(((b >> i) & 1) << perm[i] for i in range(m))
                return px | (pb << m)
            Lev.append([pos[act(g, s)] for s in X])
        # elementary transvection in GL: x_0 += x_1 on x, dual on b (b_1 += b_0)
        def g2(u):
            x = u & ((1 << m) - 1); b = u >> m
            x ^= ((x >> 1) & 1); b ^= ((b & 1) << 1)
            return x | (b << m)
        Lev.append([pos[act(g2, s)] for s in X])
        PP = Pg + Lev
        orb = -np.ones(n, dtype=int); no2 = 0
        for i in range(n):
            if orb[i] >= 0: continue
            st = [i]; orb[i] = no2
            while st:
                a = st.pop()
                for p in PP:
                    if orb[p[a]] < 0: orb[p[a]] = no2; st.append(p[a])
            no2 += 1
        bestR = (0, None)
        for mask in range(1, (1 << no2) - 1):
            U = set(i for i in range(n) if (mask >> orb[i]) & 1)
            ft = sum(1 for i in U if pt[i] not in U); fK = sum(1 for i in U if pK[i] not in U)
            if ft and fK / ft > bestR[0]: bestR = (fK / ft, (mask, fK, ft))
        print("   Theorem R, unions of %d orbits of %s: best f(K)/f(t) = %.4f %s" % (no2, pname, bestR[0], bestR[1]))
    sys.stdout.flush()
