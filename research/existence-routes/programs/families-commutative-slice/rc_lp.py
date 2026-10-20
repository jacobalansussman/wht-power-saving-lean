"""rc_lp.py (ex2-families, session 2): the slice LP strengthened by rank-cover inequalities (PROOF.md section 17).

Normalised to W = 1 (everything is homogeneous of degree 1 in (c, W)). Variables: c_z >= 0 (moves per array in
direction z). Constraints:
 (D) deletion: for every odd form Q, c(D_Q) >= 1.
 (RC) for every even form F whose support is a 2-flat (block 4) and every partition pi of its 4 points into blocks of
      size <= 2: sum over blocks B of y_B + sum over unordered cross pairs {x, y} of 2 y^C_xy >= c(F), where
      y_B <= c(B), y_B <= c(G) - 1 for every odd G whose support contains B (rank of a principal piece);
      y^C_xy <= c_x, <= c_y, <= c(G) - 1 for every odd G containing x and y (rank of a cross cell).
   Valid because E[D_F, D_F] is nonsingular of order c(F) W (Theorem 10(b)) and is the sum of its pieces, each of rank
   <= the slack of any odd support containing it (Theorem 10(c)) and <= its own size.
LP value v: every commutative-slice unit at block 4 has M >= v W. Plain is 4; the deletion LP alone gives 15/4.
usage: rc_lp.py [nopair]
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix
from refl import quad_forms

m = 4
vals, tr = quad_forms(m)
n = (1 << m) - 1
V = vals[:, 1:]                       # forms x directions (direction z -> column z-1)
odd = np.nonzero(tr == 1)[0]
even = np.nonzero(tr == 0)[0]
wt = V.sum(1)
evenF = [f for f in even if wt[f] == 4]
# check these are flats: 4 points summing to 0
evenF = [f for f in evenF if np.bitwise_xor.reduce(np.nonzero(V[f])[0] + 1) == 0]
use_pairs = not (len(sys.argv) > 1 and sys.argv[1] == "nopair")

var = {}          # name -> index
def vid(name):
    if name not in var: var[name] = len(var)
    return var[name]
for z in range(n): vid(("c", z))
rows = []         # list of (dict idx->coef, rhs) meaning sum coef x <= rhs
for q in odd:     # -c(D_Q) <= -1
    rows.append(({vid(("c", z)): -1.0 for z in np.nonzero(V[q])[0]}, -1.0))
oddsupp = [set(np.nonzero(V[q])[0].tolist()) for q in odd]
def cap_rows(yname, U, sizes):
    """y <= c(G) - 1 for odd G containing U;  y <= each size expression (list of direction sets)."""
    y = vid(yname)
    for S in oddsupp:
        if U <= S:
            d = {y: 1.0}
            for z in S: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) - 1.0
            rows.append((d, -1.0))
    for T in sizes:
        d = {y: 1.0}
        for z in T: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) - 1.0
        rows.append((d, 0.0))
    return y
for f in evenF:
    P = sorted(np.nonzero(V[f])[0].tolist())
    yS = {x: cap_rows(("S", f, x), {x}, [[x]]) for x in P}
    yP = {}
    if use_pairs:
        for a, b in itertools.combinations(P, 2):
            yP[(a, b)] = cap_rows(("P", f, a, b), {a, b}, [[a, b]])
    yC = {(a, b): cap_rows(("C", f, a, b), {a, b}, [[a], [b]]) for a, b in itertools.combinations(P, 2)}
    # partitions of P into blocks of size <= 2
    parts = [[(x,) for x in P]]
    if use_pairs:
        for a, b in itertools.combinations(P, 2):
            rest = [x for x in P if x not in (a, b)]
            parts.append([(a, b)] + [(x,) for x in rest])
        for a, b in [(P[0], P[1]), (P[0], P[2]), (P[0], P[3])]:
            rest = tuple(x for x in P if x not in (a, b))
            parts.append([(a, b), rest])
    for pi in parts:
        # c(F) - sum pieces <= 0
        d = {}
        for z in P: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) + 1.0
        blk = {}
        for B in pi:
            for x in B: blk[x] = B
            y = yS[B[0]] if len(B) == 1 else yP[tuple(sorted(B))]
            d[y] = d.get(y, 0.0) - 1.0
        for a, b in itertools.combinations(P, 2):
            if blk[a] != blk[b]:
                d[yC[(a, b)]] = d.get(yC[(a, b)], 0.0) - 2.0
        rows.append((d, 0.0))
nv = len(var)
A = lil_matrix((len(rows), nv))
bvec = np.zeros(len(rows))
for i, (d, r) in enumerate(rows):
    for j, cft in d.items(): A[i, j] = cft
    bvec[i] = r
cost = np.zeros(nv)
for z in range(n): cost[var[("c", z)]] = 1.0
bounds = [(0, None)] * n + [(None, None)] * (nv - n)
res = linprog(cost, A_ub=csr_matrix(A), b_ub=bvec, bounds=bounds, method="highs")
print(f"block 4: even 2-flats used {len(evenF)}; variables {nv}; constraints {len(rows)}; pairs {'on' if use_pairs else 'off'}")
print(f"LP status {res.status}: value {res.fun:.6f} moves per array (deletion LP alone 3.75; plain 4)")
if res.status == 0:
    c = res.x[:n]
    print("optimal c by direction:", {z + 1: round(float(c[z]), 4) for z in range(n) if c[z] > 1e-9})
