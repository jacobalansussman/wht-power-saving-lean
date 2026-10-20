#!/usr/bin/env python3
"""t_alt.py (rm2-phase): the phase-extended per-array targets at even blocks, exact checks at m = 2, 4 (and 6 for (b)-(c)).
With phases free, an array may be owed x(S) for ANY nonsingular symmetric S (two-sided free factors: x(S) lies in the big
cell P w P, as K does).  New against the no-phase game: the ALTERNATING nonsingular S (zero diagonal; even m only).
 (a) count of nonsingular alternating S (= |GL(m,2)|/|Sp(m,2)|) and of all nonsingular symmetric S, m = 2, 4;
 (b) the least k with S = sum of k matrices z z^T (rank depth, so a lone array owed x(S) needs k moves): m for S ~ I,
     m + 1 for alternating S (the z's must sum to 0);  checked by search at m = 2, 4;
 (c) deletion with all arrays owed one alternating S: on the forms Q = (a.z)(b.z) with a^T S b = 1 every array is odd and
     |supp Q| = 2^(m-2), and the orbit average gives M >= W (2^m - 1)/2^(m-2), exactly the family's cost per array.
     Checked: <Q,S> = a^T S b and |supp| = 2^(m-2) for all such pairs at m = 2, 4, 6; each z covered equally."""
import itertools
def sym_from(m, bits):
    S = [[0]*m for _ in range(m)]; k = 0
    for i in range(m):
        for j in range(i, m):
            S[i][j] = S[j][i] = (bits >> k) & 1; k += 1
    return S
def det2(S):
    m = len(S); A = [row[:] for row in S]; r = 0
    for c in range(m):
        p = next((k for k in range(r, m) if A[k][c]), None)
        if p is None: return 0
        A[r], A[p] = A[p], A[r]
        for k in range(m):
            if k != r and A[k][c]: A[k] = [x ^ y for x, y in zip(A[k], A[r])]
        r += 1
    return 1
def pair(S, a, b, m): return sum(((a >> i) & 1) * S[i][j] * ((b >> j) & 1) for i in range(m) for j in range(m)) & 1
for m in (2, 4):
    N = m * (m + 1) // 2
    allS = [sym_from(m, b) for b in range(1 << N)]
    nons = [S for S in allS if det2(S)]
    alt = [S for S in nons if all(S[i][i] == 0 for i in range(m))]
    print(f"m = {m}: nonsingular symmetric {len(nons)}, of them alternating {len(alt)}")
    # (b) rank depth by BFS over sums of z z^T
    def key(S): return tuple(tuple(r) for r in S)
    rk1 = []
    for z in range(1, 1 << m):
        v = [(z >> i) & 1 for i in range(m)]
        rk1.append(key([[v[i] & v[j] for j in range(m)] for i in range(m)]))
    depth = {key([[0]*m for _ in range(m)]): 0}; frontier = list(depth)
    k = 0
    while frontier:
        k += 1; new = []
        for S in frontier:
            for R in rk1:
                T = tuple(tuple(x ^ y for x, y in zip(r1, r2)) for r1, r2 in zip(S, R))
                if T not in depth: depth[T] = k; new.append(T)
        frontier = new
    dI = depth[key([[int(i == j) for j in range(m)] for i in range(m)])]
    dalt = sorted({depth[key(S)] for S in alt}); dnon = sorted({depth[key(S)] for S in nons if S not in alt})
    print(f"   rank depth: I -> {dI}; non-alternating nonsingular -> {dnon}; alternating nonsingular -> {dalt}")
for m in (2, 4, 6):
    S = [[0]*m for _ in range(m)]
    for i in range(0, m, 2): S[i][i+1] = S[i+1][i] = 1      # standard symplectic form
    ok = True; cover = [0] * (1 << m); npairs = 0
    for a in range(1, 1 << m):
        for b in range(1, 1 << m):
            if a == b or not pair(S, a, b, m): continue
            npairs += 1
            # Q(z) = (a.z)(b.z); its coefficients q_ij = a_i b_j + a_j b_i (i<j), q_i = a_i b_i; <Q,S> = sum_{i<j} q_ij S_ij
            q = sum((((a >> i) & 1) * ((b >> j) & 1) + ((a >> j) & 1) * ((b >> i) & 1)) * S[i][j]
                    for i in range(m) for j in range(i + 1, m)) & 1
            supp = [z for z in range(1, 1 << m) if (bin(a & z).count("1") & 1) and (bin(b & z).count("1") & 1)]
            ok &= (q == 1) and len(supp) == 2 ** (m - 2)
            for z in supp: cover[z] += 1
    eq = len(set(cover[1:])) == 1
    print(f"m = {m}: {npairs} ordered pairs (a,b) with a^T S b = 1: every Q = (a.z)(b.z) odd at S with support 2^(m-2): {ok}; "
          f"every direction covered equally ({cover[1]} times): {eq}; deletion: M >= W * (2^m-1)/2^(m-2) = {(2**m-1)/2**(m-2)} W")
