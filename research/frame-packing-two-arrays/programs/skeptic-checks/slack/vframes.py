#!/usr/bin/env python3
"""v-slack own check of the frame facts both claims use, written from scratch (no import of their code).

(1) From actual matrices (exact Gaussian rationals) at m = 1, 2, 3: the conjugation action of the engine letters
    on Pauli labels (a|b).  dir z = (1+i)/2 I + (1-i)/2 X^z;  free letters: X^z, Z^z, address maps (GL(m,2)
    permutations of the basis index), phase letters diag(i^{(z.x) mod 2}).
(2) Lagrangians of F_2^{2m} (dual polar graph), L0 = Z-type {(0|b)}: free letters fix L0, each dir z moves every
    Lagrangian by distance <= 1 and L0 by exactly 1, d(L0, Kbar L0) = m with K = dir e1 ... dir em.
(3) Triangles / lateral steps exist (so path-length parity is NOT fixed: relevant to the parity question).
(4) Cut constant c_m = max_S [share of opposite pairs separated] / [share of edges cut]:
    brute force over all subsets at m = 2; point cuts {L : v in L} at m = 3.
"""
from fractions import Fraction as F
from itertools import product, combinations
import json, sys

# ---------- exact Gaussian rationals as (re, im) ----------
def cadd(x, y): return (x[0] + y[0], x[1] + y[1])
def cmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
ZERO = (F(0), F(0)); ONE = (F(1), F(0)); I = (F(0), F(1))

def matmul(A, B):
    n = len(A); p = len(B[0]); k = len(B)
    C = [[ZERO] * p for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for t in range(k):
            a = Ai[t]
            if a == ZERO: continue
            Bt = B[t]
            for j in range(p):
                if Bt[j] != ZERO:
                    C[i][j] = cadd(C[i][j], cmul(a, Bt[j]))
    return C

def conjT(A):
    return [[(A[j][i][0], -A[j][i][1]) for j in range(len(A))] for i in range(len(A[0]))]

def ident(n): return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]

def dot(a, b): return bin(a & b).count('1') & 1

def Xz(m, z):
    n = 1 << m
    return [[ONE if i == (j ^ z) else ZERO for j in range(n)] for i in range(n)]

def Zz(m, z):
    n = 1 << m
    return [[((F((-1) ** dot(z, i)), F(0)) if i == j else ZERO) for j in range(n)] for i in range(n)]

def pauli(m, a, b):  # X^a Z^b
    return matmul(Xz(m, a), Zz(m, b))

def dirz(m, z):
    ca = (F(1, 2), F(1, 2)); cb = (F(1, 2), F(-1, 2))
    n = 1 << m; X = Xz(m, z)
    return [[cadd(ca if i == j else ZERO, cmul(cb, X[i][j])) for j in range(n)] for i in range(n)]

def phase(m, z):  # diag(i^{(z.x) mod 2})
    n = 1 << m
    return [[(I if dot(z, i) else ONE) if i == j else ZERO for j in range(n)] for i in range(n)]

def addrmap(m, Amat):  # permutation |x> -> |A x>, A an m x m invertible bit matrix given as list of column ints
    n = 1 << m
    def ap(x):
        y = 0
        for c in range(m):
            if (x >> c) & 1: y ^= Amat[c]
        return y
    return [[ONE if i == ap(j) else ZERO for j in range(n)] for i in range(n)]

def ident_pauli(m, M):
    """M must be a phase times a Pauli X^a Z^b; return (a, b)."""
    n = 1 << m
    for a in range(n):
        for b in range(n):
            P = pauli(m, a, b)
            # find ratio
            ratio = None; ok = True
            for i in range(n):
                for j in range(n):
                    if (P[i][j] == ZERO) != (M[i][j] == ZERO): ok = False; break
                    if P[i][j] != ZERO:
                        # M = c P -> c = M/P, P entries are +-1
                        c = (M[i][j][0] * P[i][j][0], M[i][j][1] * P[i][j][0])
                        if ratio is None: ratio = c
                        elif ratio != c: ok = False; break
                if not ok: break
            if ok and ratio is not None: return (a, b)
    raise ValueError("not a Pauli multiple")

def sp_action(m, g):
    """images of the basis Paulis (e_k|0), (0|e_k) under P -> g P g^-1; returns a function on (a|b) ints."""
    ginv = conjT(g)  # all our letters are unitary
    imgs = {}
    for k in range(m):
        for (a, b) in [((1 << k), 0), (0, (1 << k))]:
            imgs[(a, b)] = ident_pauli(m, matmul(matmul(g, pauli(m, a, b)), ginv))
    def act(v):
        a, b = v; ra, rb = 0, 0
        for k in range(m):
            if (a >> k) & 1: x = imgs[((1 << k), 0)]; ra ^= x[0]; rb ^= x[1]
            if (b >> k) & 1: x = imgs[(0, (1 << k))]; ra ^= x[0]; rb ^= x[1]
        return (ra, rb)
    return act

def symp(u, v): return dot(u[0], v[1]) ^ dot(u[1], v[0])

def span(vecs):
    S = {(0, 0)}
    for v in vecs:
        S |= {(s[0] ^ v[0], s[1] ^ v[1]) for s in S}
    return frozenset(S)

def lagrangians(m):
    """all m-dim totally isotropic subspaces of F_2^{2m} (brute force over spans of isotropic sets via BFS)."""
    n = 1 << m
    allv = [(a, b) for a in range(n) for b in range(n) if (a, b) != (0, 0)]
    found = set(); L0 = span([(0, 1 << k) for k in range(m)])
    # BFS under all transvections (they generate Sp(2m,2), which is transitive on Lagrangians)
    def tv(v, L):
        return frozenset(((x[0] ^ v[0], x[1] ^ v[1]) if symp(x, v) else x) for x in L)
    frontier = [L0]; found.add(L0)
    while frontier:
        nf = []
        for L in frontier:
            for v in allv:
                L2 = tv(v, L)
                if L2 not in found: found.add(L2); nf.append(L2)
        frontier = nf
    return L0, sorted(found, key=lambda s: sorted(s))

def dist(m, A, B):
    k = len(A & B).bit_length() - 1
    return m - k

def check_m(m, out):
    n = 1 << m
    res = {}
    L0, Ls = lagrangians(m)
    res['n_lagrangians'] = len(Ls)
    # isotropy sanity
    res['all_isotropic'] = all(symp(u, v) == 0 for L in Ls for u in L for v in L)
    idx = {L: i for i, L in enumerate(Ls)}
    # letters from actual matrices
    letters_free = {}
    for z in range(1, n):
        letters_free['X^%d' % z] = Xz(m, z)
        letters_free['Z^%d' % z] = Zz(m, z)
        letters_free['ph%d' % z] = phase(m, z)
    # some address maps: all elementary transvection matrices x_i += x_j and swaps
    for i in range(m):
        for j in range(m):
            if i != j:
                cols = [1 << c for c in range(m)]; cols[j] ^= (1 << i)   # e_j -> e_j + e_i
                letters_free['add%d%d' % (i, j)] = addrmap(m, cols)
    paid = {('dir%d' % z): dirz(m, z) for z in range(1, n)}
    acts_free = {k: sp_action(m, g) for k, g in letters_free.items()}
    acts_paid = {k: sp_action(m, g) for k, g in paid.items()}
    def img(act, L): return frozenset(act(v) for v in L)
    res['free_fix_L0'] = all(img(a, L0) == L0 for a in acts_free.values())
    res['paid_move_L0_exactly_1'] = all(dist(m, img(a, L0), L0) == 1 for a in acts_paid.values())
    res['paid_move_every_L_at_most_1'] = all(dist(m, img(a, L), L) <= 1 for a in acts_paid.values() for L in Ls)
    # the transvection formula (a|b) -> (a + (z.b) z | b) for dir z
    ok = True
    for z in range(1, n):
        a_ = acts_paid['dir%d' % z]
        for a in range(n):
            for b in range(n):
                if a_((a, b)) != (a ^ (z if dot(z, b) else 0), b): ok = False
    res['dir_is_transvection_(z|0)'] = ok
    # K = dir e1 ... dir em
    K = ident(n)
    for k in range(m): K = matmul(K, dirz(m, 1 << k))
    aK = sp_action(m, K)
    res['d(L0, Kbar L0)'] = dist(m, img(aK, L0), L0)
    # graph facts
    N = len(Ls)
    D = [[dist(m, A, B) for B in Ls] for A in Ls]
    edges = [(i, j) for i in range(N) for j in range(i + 1, N) if D[i][j] == 1]
    opp = [(i, j) for i in range(N) for j in range(i + 1, N) if D[i][j] == m]
    res['n_edges'] = len(edges); res['n_opposite_pairs'] = len(opp)
    # triangles: a_1 = number of common neighbours of an adjacent pair
    i, j = edges[0]
    res['a1_common_neighbours_of_an_edge'] = sum(1 for k in range(N) if D[i][k] == 1 and D[j][k] == 1)
    # cut constant
    nb = {v: set() for v in range(N)}
    for (u, v) in edges: nb[u].add(v); nb[v].add(u)
    if m == 2:
        best = F(0); bestS = None; cnt = 0
        for mask in range(1, 1 << (N - 1)):  # vertex N-1 always outside: subsets up to complement
            cnt += 1
            S = [(mask >> v) & 1 for v in range(N)]
            ce = sum(1 for (u, v) in edges if S[u] != S[v])
            co = sum(1 for (u, v) in opp if S[u] != S[v])
            if ce == 0: continue
            r = F(co, len(opp)) / F(ce, len(edges))
            if r > best: best = r; bestS = [v for v in range(N) if S[v]]
        res['c2_bruteforce_max_ratio'] = str(best); res['c2_subsets_checked'] = cnt
        res['c2_best_set_size'] = len(bestS)
    # point cuts {L : v in L}
    pts = []
    for a in range(n):
        for b in range(n):
            if (a, b) == (0, 0): continue
            S = [1 if (a, b) in L else 0 for L in Ls]
            ce = sum(1 for (u, v) in edges if S[u] != S[v])
            co = sum(1 for (u, v) in opp if S[u] != S[v])
            pts.append(F(co, len(opp)) / F(ce, len(edges)))
    res['point_cut_ratios'] = sorted(set(str(x) for x in pts))
    res['formula_2(1-2^-m)'] = str(2 * (1 - F(1, 2 ** m)))
    out[m] = res
    return res

if __name__ == '__main__':
    out = {}
    for m in (1, 2, 3):
        r = check_m(m, out)
        print('m =', m)
        for k, v in r.items(): print('   ', k, ':', v)
        sys.stdout.flush()
    json.dump({str(k): v for k, v in out.items()}, open(__file__.replace('.py', '.json'), 'w'), indent=1)
