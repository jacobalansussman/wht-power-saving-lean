"""exact_cert.py: turn a floating HiGHS optimum of  min cost^T x, A x <= b, x_j >= 0 (j < n), x_j free (j >= n)
into an EXACT rational dual certificate lambda >= 0 with A_free^T lambda = 0, cost_c + A_c^T lambda >= 0, and report
-b^T lambda (a proved lower bound on the primal value).  A must have small integer entries."""
from fractions import Fraction as Fr
import numpy as np

def certify(A, b, cost, n, res, target, tol=1e-9):
    nv = A.shape[1]
    lam = -res.ineqlin.marginals
    sup = [i for i in range(A.shape[0]) if lam[i] > tol]
    As = A[sup].toarray(); bs = b[sup]
    act_c = [j for j in range(n) if abs(cost[j] + As[:, j] @ lam[sup]) < 1e-7]
    touched = [j for j in range(n, nv) if np.any(As[:, j] != 0)]
    Mx = [[Fr(int(round(As[i][j]))) for i in range(len(sup))] for j in touched]
    rhs = [Fr(0)] * len(touched)
    for j in act_c:
        Mx.append([Fr(int(round(As[i][j]))) for i in range(len(sup))]); rhs.append(Fr(-int(round(cost[j]))))
    Mx.append([Fr(-int(round(bs[i]))) for i in range(len(sup))]); rhs.append(Fr(target))
    nrow, ncol = len(Mx), len(sup)
    aug = [row[:] + [rhs[k]] for k, row in enumerate(Mx)]
    piv = []; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, nrow) if aug[i][c] != 0), None)
        if p is None: continue
        aug[r], aug[p] = aug[p], aug[r]
        inv = 1 / aug[r][c]; aug[r] = [v * inv for v in aug[r]]
        for i in range(nrow):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]; aug[i] = [vi - f * vr for vi, vr in zip(aug[i], aug[r])]
        piv.append(c); r += 1
        if r == nrow: break
    incons = any(all(v == 0 for v in aug[i][:ncol]) and aug[i][ncol] != 0 for i in range(r, nrow))
    free = [c for c in range(ncol) if c not in piv]
    sol = [Fr(0)] * ncol
    for c in free: sol[c] = Fr(float(lam[sup[c]])).limit_denominator(10**6)
    for k, c in enumerate(piv): sol[c] = aug[k][ncol] - sum(aug[k][j] * sol[j] for j in free)
    # exact verification over all columns
    Ac = A[sup].tocoo()
    colsum = {}
    for i, j, v in zip(Ac.row, Ac.col, Ac.data):
        colsum[j] = colsum.get(j, Fr(0)) + sol[i] * Fr(int(round(v)))
    ok_nonneg = all(v >= 0 for v in sol)
    ok_free = all(colsum.get(j, 0) == 0 for j in range(n, nv))
    ok_c = all(Fr(int(round(cost[j]))) + colsum.get(j, 0) >= 0 for j in range(n))
    val = -sum(sol[k] * Fr(int(round(bs[k]))) for k in range(ncol))
    good = ok_nonneg and ok_free and ok_c and not incons
    return dict(support=len(sup), free_vars=len(free), inconsistent=incons, nonneg=ok_nonneg, free_cols_cancel=ok_free,
                c_cols_ok=ok_c, value=val, certificate=good, sol=sol, sup=sup)
