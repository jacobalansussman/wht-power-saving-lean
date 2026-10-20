"""Exact checks behind the dense-pair lemma (07-pivots.tex) and its extension.

Run: python3 -I dense_pair.py   (pure Python, exact rationals, a few seconds)

1. Re-derive the paper's 4x4 identity: L2 L1 T R1 R2 equals B(u)(x)B(u) in the stated row/column
   orders, with T = I + Delta k^T/(1-lambda), u = 1/lambda.  [calibration]
2. Ask whether B(u)^(x)3 (8x8) is, like the 4x4 case, a permutation-and-scaling of
   "diagonal + rank one" (the property that lets four exchanges and two shears reach it from the
   4-cycle I+J). Entries are u^(alpha.gamma); for generic u, off-diagonal rank one after scaling
   means exponent(i,j) = a_i + b_j on every off-diagonal position, i.e. for every pair of columns
   the exponent difference is the same on every row off both their diagonal positions.
   Only the column permutation matters (simultaneous permutations keep the diagonal).
"""
from fractions import Fraction as Fr
from itertools import permutations, product

# ---- 1. the 4x4 identity -------------------------------------------------------------------
lam = Fr(3, 7)
u = 1 / lam
k = [1, -1, 1, -1]
Delta = [-1, 1, lam, -lam]
T = [[(1 if i == j else 0) + Delta[i] * k[j] / (1 - lam) for j in range(4)] for i in range(4)]
L1 = [(1 - lam) / Delta[i] for i in range(4)]
R1 = [1, -1, 1, -1]
L2 = [u, u, u, 1]
R2 = [1, 1, 1, lam]
S = [[L2[i] * L1[i] * T[i][j] * R1[j] * R2[j] for j in range(4)] for i in range(4)]
rows = [(1, 0), (0, 1), (1, 1), (0, 0)]
cols = [(0, 1), (1, 0), (1, 1), (0, 0)]
BB = [[u ** (a[0] * g[0] + a[1] * g[1]) for g in cols] for a in rows]
print("4x4 identity L2 L1 T R1 R2 == B(u)(x)B(u):", S == BB)

# ---- 2. diagonal + rank one test ------------------------------------------------------------
def exps(nbits):
    pts = list(product((0, 1), repeat=nbits))
    return [[sum(x * y for x, y in zip(a, g)) for g in pts] for a in pts]

def has_diag_plus_rank_one(E):
    n = len(E)
    hits = 0
    for perm in permutations(range(n)):
        ok = True
        for j in range(n):
            if not ok:
                break
            for l in range(j + 1, n):
                diffs = {E[i][perm[j]] - E[i][perm[l]] for i in range(n) if i != j and i != l}
                if len(diffs) > 1:
                    ok = False
                    break
        if ok:
            hits += 1
    return hits

print("2 bits: column orders giving diagonal + rank one:", has_diag_plus_rank_one(exps(2)), "of 24")
print("3 bits: column orders giving diagonal + rank one:", has_diag_plus_rank_one(exps(3)), "of 40320")
