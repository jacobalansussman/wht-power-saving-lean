"""rc_lp2.py: the block-4 rank-cover LP in its smallest form (no pair pieces), with an EXACT rational certificate.

Primal (W = 1): min sum_z c_z over c >= 0 and free y's:
  (D_Q)   -c(D_Q) <= -1                           for odd Q
  (S_Fx)  y_Fx - c_x <= 0 ; y_Fx - c(D_G) <= -1   for odd G with x in D_G        (diagonal block of E_F at x)
  (C_Fxy) z_Fxy - c_x <= 0, z_Fxy - c_y <= 0 ; z_Fxy - c(D_G) <= -1 for odd G containing x, y   (cross block, both orders)
  (R_F)   c(F) - sum_x y_Fx - 2 sum_{x<y} z_Fxy <= 0                             for even 2-flats F
Step 1: floating LP (HiGHS) and its dual. Step 2: keep the rows with dual > 1e-9, re-solve the dual restricted to them
EXACTLY (fractions, by solving the tight equality system), and verify in exact arithmetic that lambda >= 0, A_y^T lambda = 0,
cost_c + A_c^T lambda >= 0 and -b^T lambda = 4.  That proves M >= 4W for every commutative-slice unit at block 4.
"""
import sys, itertools, json
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from refl import quad_forms

m = 4
vals, tr = quad_forms(m)
n = (1 << m) - 1
V = vals[:, 1:]
odd = [int(q) for q in np.nonzero(tr == 1)[0]]
wt = V.sum(1)
evenF = [int(f) for f in np.nonzero(tr == 0)[0] if wt[f] == 4 and np.bitwise_xor.reduce(np.nonzero(V[f])[0] + 1) == 0]
supp = {q: frozenset(np.nonzero(V[q])[0].tolist()) for q in range(len(V))}
oddsupp = sorted(set(supp[q] for q in odd), key=lambda s: (len(s), sorted(s)))

var = {}
def vid(k):
    if k not in var: var[k] = len(var)
    return var[k]
for z in range(n): vid(("c", z))
rows, labels = [], []
def add(d, r, lab): rows.append((d, r)); labels.append(lab)
for S in oddsupp:
    add({vid(("c", z)): -1 for z in S}, -1, ("D", sorted(S)))
for f in evenF:
    P = sorted(supp[f])
    for x in P:
        y = vid(("y", f, x))
        add({y: 1, vid(("c", x)): -1}, 0, ("Ssize", f, x))
        for S in oddsupp:
            if x in S:
                d = {y: 1}
                for z in S: d[vid(("c", z))] = d.get(vid(("c", z)), 0) - 1
                add(d, -1, ("Sslack", f, x, sorted(S)))
    for a, b in itertools.combinations(P, 2):
        zz = vid(("z", f, a, b))
        add({zz: 1, vid(("c", a)): -1}, 0, ("Csize", f, a, b, a))
        add({zz: 1, vid(("c", b)): -1}, 0, ("Csize", f, a, b, b))
        for S in oddsupp:
            if a in S and b in S:
                d = {zz: 1}
                for z in S: d[vid(("c", z))] = d.get(vid(("c", z)), 0) - 1
                add(d, -1, ("Cslack", f, a, b, sorted(S)))
    d = {}
    for x in P: d[vid(("c", x))] = d.get(vid(("c", x)), 0) + 1
    for x in P: d[vid(("y", f, x))] = -1
    for a, b in itertools.combinations(P, 2): d[vid(("z", f, a, b))] = -2
    add(d, 0, ("R", f))
nv = len(var); nr = len(rows)
ri, ci, vv = [], [], []
for i, (d, r) in enumerate(rows):
    for j, cft in d.items(): ri.append(i); ci.append(j); vv.append(cft)
A = csr_matrix((vv, (ri, ci)), shape=(nr, nv)); b = np.array([r for d, r in rows], float)
cost = np.zeros(nv); cost[:n] = 1
res = linprog(cost, A_ub=A, b_ub=b, bounds=[(0, None)] * n + [(None, None)] * (nv - n), method="highs")
lam = -res.ineqlin.marginals       # HiGHS: marginals <= 0 for <= rows
print(f"variables {nv}, rows {nr}; LP value {res.fun:.9f}; dual objective {-(b @ lam):.9f}; dual support {int((lam > 1e-9).sum())}")
sup = [i for i in range(nr) if lam[i] > 1e-9]
# exact: the restricted dual  max -b^T l  s.t. A_y^T l = 0, cost_c + A_c^T l >= 0, l >= 0  on the support rows.
# Solve it exactly with a rational LP via fractions: use the floating optimum's active c-columns as equalities.
ycols = list(range(n, nv))
As = A[sup].toarray()
act_c = [j for j in range(n) if abs(cost[j] + As[:, j] @ lam[sup]) < 1e-9]
# equality system: A_y^T l = 0 (all y columns touched by sup), cost_c + A_c^T l = 0 for active c, -b^T l = 4
touched = [j for j in ycols if np.any(As[:, j] != 0)]
M = [[Fr(int(round(As[i][j]))) for i in range(len(sup))] for j in touched]
rhs = [Fr(0)] * len(touched)
for j in act_c:
    M.append([Fr(int(round(As[i][j]))) for i in range(len(sup))]); rhs.append(Fr(-1))
M.append([Fr(-int(round(b[sup[i]]))) for i in range(len(sup))]); rhs.append(Fr(4))
# least-norm-free exact solve: Gaussian elimination, free variables set from the floating solution (rationalised)
import copy
nrow, ncol = len(M), len(sup)
aug = [row[:] + [rhs[k]] for k, row in enumerate(M)]
piv = []; r = 0
for c in range(ncol):
    p = next((i for i in range(r, nrow) if aug[i][c] != 0), None)
    if p is None: continue
    aug[r], aug[p] = aug[p], aug[r]
    inv = 1 / aug[r][c]
    aug[r] = [v * inv for v in aug[r]]
    for i in range(nrow):
        if i != r and aug[i][c] != 0:
            fct = aug[i][c]; aug[i] = [vi - fct * vr for vi, vr in zip(aug[i], aug[r])]
    piv.append(c); r += 1
    if r == nrow: break
incons = any(all(v == 0 for v in aug[i][:ncol]) and aug[i][ncol] != 0 for i in range(r, nrow))
free = [c for c in range(ncol) if c not in piv]
lf = {c: Fr(lam[sup[c]]).limit_denominator(64) for c in free}
sol = [Fr(0)] * ncol
for c in free: sol[c] = lf[c]
for k, c in enumerate(piv):
    sol[c] = aug[k][ncol] - sum(aug[k][j] * sol[j] for j in free)
# exact verification on the FULL row set (lambda = 0 off the support)
L = {sup[k]: sol[k] for k in range(ncol)}
ok_nonneg = all(v >= 0 for v in L.values())
colsum = {}
for i, l in L.items():
    for j, cft in rows[i][0].items():
        colsum[j] = colsum.get(j, Fr(0)) + l * cft
ok_y = all(colsum.get(j, 0) == 0 for j in ycols)
ok_c = all(Fr(1) + colsum.get(j, 0) >= 0 for j in range(n))
obj = -sum(l * rows[i][1] for i, l in L.items())
print(f"exact certificate: inconsistent={incons}; lambda >= 0: {ok_nonneg}; y-columns cancel: {ok_y}; c-columns ok: {ok_c}; "
      f"value {obj}  -> {'PROVES M >= 4W' if (ok_nonneg and ok_y and ok_c and obj == 4 and not incons) else 'not a certificate'}")
kinds = {}
for i, l in L.items():
    if l != 0: kinds[labels[i][0]] = kinds.get(labels[i][0], 0) + 1
print("rows used by kind:", kinds)
out = [{"row": labels[i], "lambda": str(l)} for i, l in L.items() if l != 0]
json.dump(out, open("rc_cert_b4.json", "w"), default=str, indent=0)
print("certificate written to rc_cert_b4.json")
