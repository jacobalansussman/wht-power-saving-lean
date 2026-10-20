"""f4_label.py (g-small, H1/H2/H13): decoders B = U A U^T in the LABEL family over F_2 and F_4.
B[S,T] = u_S^T A u_T must vanish on non-orthogonal E8 pairs; diag u_S^T A u_S != 0.  The constraint is F_2-linear in A,
so the F_4 solutions are W (x) F_4 with W the F_2 solution space.  A scalar F_4 decoder of rank 8 would be an F_2
vector-linear (t = 2, two calls sharing frames) decoder with 16 terms per 2 calls: cst 64 per call instead of 72."""
import sys, itertools; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
B0 = [[(1 + G[i][j]) % 2 for j in range(V)] for i in range(V)]
U = cf.f2_labels(B0)
# unknowns A[a][b], 81 of them; constraint rows: for non-orth (S,T): sum_ab u_S[a] u_T[b] A[a][b] = 0
rows = set()
for S in range(V):
    for T in range(V):
        if NONORTH[S][T]:
            rows.add(sum(1 << (9 * a + b) for a in range(9) for b in range(9) if U[S][a] and U[T][b]))
rows = list(rows)
# nullspace over F_2 (bitmask gaussian elimination)
piv = {}
for r in rows:
    x = r
    while x:
        h = x.bit_length() - 1
        if h in piv: x ^= piv[h]
        else: piv[h] = x; break
# reduce fully
keys = sorted(piv)
for h in keys:
    for h2 in keys:
        if h2 != h and piv[h2] >> h & 1: piv[h2] ^= piv[h]
free = [c for c in range(81) if c not in piv]
basis = []
for f in free:
    x = 1 << f
    for h, r in piv.items():
        if r >> f & 1: x |= 1 << h
    basis.append(x)
print('constraint rank', len(piv), ' dim W (F_2 solutions A) =', len(basis))
def mat(x): return [[x >> (9 * a + b) & 1 for b in range(9)] for a in range(9)]
for x in basis: 
    M = mat(x); print('  basis element: rank_F2', cf.rank(M, 2), ' symmetric', all(M[a][b] == M[b][a] for a in range(9) for b in range(9)),
          ' diag-ok on all lines', all(sum(U[S][a] * M[a][b] * U[S][b] for a in range(9) for b in range(9)) % 2 for S in range(V)))
# F_4 = F_2[w]/(w^2+w+1): elements 0,1,w,w^2 = (c0,c1) as c0 + c1 w.  Enumerate A = sum_i lam_i W_i, lam in F_4.
def f4mul(a, b):  # a=(a0,a1)
    a0, a1 = a; b0, b1 = b
    # (a0 + a1 w)(b0 + b1 w) = a0b0 + (a0b1+a1b0) w + a1b1 w^2, w^2 = w + 1
    return ((a0 & b0) ^ (a1 & b1), (a0 & b1) ^ (a1 & b0) ^ (a1 & b1))
def f4rank(M):
    M = [list(r) for r in M]; n = len(M); m = len(M[0]); rk = 0
    inv = {(1, 0): (1, 0), (0, 1): (1, 1), (1, 1): (0, 1)}
    for c in range(m):
        p = next((i for i in range(rk, n) if M[i][c] != (0, 0)), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        iv = inv[M[rk][c]]
        M[rk] = [f4mul(iv, e) for e in M[rk]]
        for i in range(n):
            if i != rk and M[i][c] != (0, 0):
                f = M[i][c]
                M[i] = [(e[0] ^ g[0], e[1] ^ g[1]) for e, g in zip(M[i], [f4mul(f, e2) for e2 in M[rk]])]
        rk += 1
    return rk
E = [(0, 0), (1, 0), (0, 1), (1, 1)]
k = len(basis); mats = [mat(x) for x in basis]
best = {}
cnt = 0
for lam in itertools.product(E, repeat=k):
    if all(l == (0, 0) for l in lam): continue
    A = [[(0, 0)] * 9 for _ in range(9)]
    for l, M in zip(lam, mats):
        for a in range(9):
            for b in range(9):
                if M[a][b]: A[a][b] = (A[a][b][0] ^ l[0], A[a][b][1] ^ l[1])
    # diag condition: u^T A u = sum_{a,b in supp u} A[a][b] != 0
    ok = True
    for S in range(V):
        s = (0, 0)
        su = [a for a in range(9) if U[S][a]]
        for a in su:
            for b in su: s = (s[0] ^ A[a][b][0], s[1] ^ A[a][b][1])
        if s == (0, 0): ok = False; break
    cnt += 1
    if ok:
        r = f4rank(A); best[r] = best.get(r, 0) + 1
print('F_4 combinations tried', cnt, ' fitting decoders by rank:', dict(sorted(best.items())))
