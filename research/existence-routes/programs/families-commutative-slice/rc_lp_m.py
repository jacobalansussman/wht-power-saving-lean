"""rc_lp_m.py: the rank-cover LP at block m (4 or 5), singleton + cross-cell pieces (Theorem 17's form), compact version.

W = 1. Variables c_z >= 0; y_x (bound for the rank of E[x-moves, x-moves]) and z_xy (bound for the ranks of E[x, y] and
E[y, x]); these bounds do not depend on the even form, so one variable each:
  y_x <= c_x,  y_x <= c(D_G) - 1 for odd G with x in D_G;
  z_xy <= c_x, z_xy <= c_y, z_xy <= c(D_G) - 1 for odd G containing x and y;
  (R_F) c(F) <= sum_{x in F} y_x + 2 sum_{x<y in F} z_xy   for every even form F with |supp F| <= maxF;
  (D_Q) c(D_Q) >= 1 for every odd Q.
Odd forms used in the bounds: those with weight <= maxG (fewer forms = weaker but still valid bounds).
usage: rc_lp_m.py m maxG maxF
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from refl import quad_forms

m, maxG, maxF = map(int, sys.argv[1:4])
vals, tr = quad_forms(m)
V = vals[:, 1:].astype(np.int8); n = (1 << m) - 1
wt = V.sum(1)
oddS = sorted(set(frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero(tr == 1)[0]), key=len)
evenS = sorted(set(frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero((tr == 0) & (wt > 0) & (wt <= maxF))[0]), key=len)
boundS = [S for S in oddS if len(S) <= maxG]
iy = {x: n + x for x in range(n)}
pairs = list(itertools.combinations(range(n), 2))
iz = {p: 2 * n + k for k, p in enumerate(pairs)}
nv = 2 * n + len(pairs)
ri, ci, vv, b = [], [], [], []
r = 0
def row(d, rhs):
    global r
    for j, cft in d.items(): ri.append(r); ci.append(j); vv.append(cft)
    b.append(rhs); r += 1
for S in oddS:                                   # deletion
    row({z: -1.0 for z in S}, -1.0)
for x in range(n):
    row({iy[x]: 1.0, x: -1.0}, 0.0)
for (a, c) in pairs:
    row({iz[(a, c)]: 1.0, a: -1.0}, 0.0); row({iz[(a, c)]: 1.0, c: -1.0}, 0.0)
for S in boundS:
    L = sorted(S)
    for x in L:
        d = {z: -1.0 for z in L}; d[iy[x]] = 1.0; row(d, -1.0)
    for (a, c) in itertools.combinations(L, 2):
        d = {z: -1.0 for z in L}; d[iz[(a, c)]] = 1.0; row(d, -1.0)
nR = 0
for F in evenS:
    L = sorted(F)
    d = {z: 1.0 for z in L}
    for x in L: d[iy[x]] = -1.0
    for (a, c) in itertools.combinations(L, 2): d[iz[(a, c)]] = -2.0
    row(d, 0.0); nR += 1
A = csr_matrix((vv, (ri, ci)), shape=(r, nv))
cost = np.zeros(nv); cost[:n] = 1
res = linprog(cost, A_ub=A, b_ub=np.array(b), bounds=[(0, None)] * n + [(None, None)] * (nv - n), method="highs")
print(f"block {m}: odd supports {len(oddS)} (bounds from weight <= {maxG}: {len(boundS)}); even supports used "
      f"(weight <= {maxF}): {nR}; rows {r}; LP status {res.status}; value {res.fun if res.status == 0 else None} "
      f"(plain {m}; deletion LP {(2**m - 1) / 2**(m - 2)})", flush=True)
if res.status == 0:
    c = res.x[:n]
    print("optimal c:", {z + 1: round(float(c[z]), 4) for z in range(n) if c[z] > 1e-9})
