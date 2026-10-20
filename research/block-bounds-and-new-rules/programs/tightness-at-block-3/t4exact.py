#!/usr/bin/env python3
"""t4exact.py (OPTIONAL item, g-rules R2): exact rank over Q(i) of the single-fibre matrix of g-rules/t4_partial_cap.py at
m=3, f=5:  B = pi (2M) (x) Lambda^4 (2M) |_N  (rows: basis e_x (x) e_T with x not in T; columns: x in T; 280 x 280),
2M = (1+i) on the diagonal of the pair {0,1}, (1-i) off-diagonal on it, 2 on the other six coordinates.
Exact arithmetic: Gaussian integers as Python int pairs; elimination over Q(i) with Fractions.  No floats, no mod p.
"""
import itertools
from fractions import Fraction
n = 8; k = 4
S4 = list(itertools.combinations(range(n), k))
basis = [(x, T) for x in range(n) for T in S4]
inN = [x in T for (x, T) in basis]
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gadd(a, b): return (a[0] + b[0], a[1] + b[1])
# 2M as a dict of columns: col c -> list of (row, gaussian int)
M2 = {c: [(c, (2, 0))] for c in range(n)}
M2[0] = [(0, (1, 1)), (1, (1, -1))]; M2[1] = [(1, (1, 1)), (0, (1, -1))]
def sort_sign(idx):
    s = 1; idx = list(idx)
    for i in range(len(idx)):
        for j in range(len(idx) - 1 - i):
            if idx[j] > idx[j + 1]: idx[j], idx[j + 1] = idx[j + 1], idx[j]; s = -s
    return s, tuple(idx)
def wedge_apply(T):
    out = {}
    for choice in itertools.product(*[M2[t] for t in T]):
        rows = [r for (r, c) in choice]
        if len(set(rows)) < k: continue
        coef = (1, 0)
        for (r, c) in choice: coef = gm(coef, c)
        s, key = sort_sign(rows)
        out[key] = gadd(out.get(key, (0, 0)), (s * coef[0], s * coef[1]))
    return {kk: v for kk, v in out.items() if v != (0, 0)}
cols = []
for j, (x, T) in enumerate(basis):
    if not inN[j]: continue
    wT = wedge_apply(T); col = {}
    for (y, cy) in M2[x]:
        for T2, cT in wT.items():
            if y in T2: continue          # keep only the rows outside N (projection pi)
            key = (y, T2); col[key] = gadd(col.get(key, (0, 0)), gm(cy, cT))
    cols.append({kk: v for kk, v in col.items() if v != (0, 0)})
nz = [c for c in cols if c]
rows = sorted({r for c in nz for r in c})
print(f"B: {sum(1 for b in basis if not (b[0] in b[1]))} x {len(cols)}; nonzero columns {len(nz)}, nonzero rows {len(rows)}")
# exact elimination over Q(i): complex numbers as (Fraction, Fraction)
def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def csub(a, b): return (a[0] - b[0], a[1] - b[1])
def cinv(a):
    d = a[0] * a[0] + a[1] * a[1]; return (a[0] / d, -a[1] / d)
ridx = {r: i for i, r in enumerate(rows)}
mat = [[(Fraction(0), Fraction(0))] * len(rows) for _ in nz]      # one row per nonzero column of B (rank of B^T)
for i, c in enumerate(nz):
    for r, v in c.items(): mat[i][ridx[r]] = (Fraction(v[0]), Fraction(v[1]))
rank = 0; ncol = len(rows); R = mat
for c in range(ncol):
    piv = next((i for i in range(rank, len(R)) if R[i][c] != (0, 0)), None)
    if piv is None: continue
    R[rank], R[piv] = R[piv], R[rank]
    inv = cinv(R[rank][c]); R[rank] = [cmul(v, inv) for v in R[rank]]
    for i in range(rank + 1, len(R)):
        if R[i][c] != (0, 0):
            fct = R[i][c]; R[i] = [csub(v, cmul(fct, w)) for v, w in zip(R[i], R[rank])]
    rank += 1
print(f"EXACT rank over Q(i) of the single-fibre matrix B (m=3, f=5): {rank}")
