# checks for refl.py: (1) F_Q formula = compression of the ordered product of reflections, (2) Jacobian vs finite
# differences, (3) the plain program in (A,B,H) form has zero residual.
import numpy as np
from refl import System, quad_forms, deletion_bound
rng = np.random.default_rng(1)
m, W = 3, 2
S = [1, 2, 4, 3, 5, 6, 7]
sy = System(m, S, W)
x = (rng.standard_normal(sy.nvar) + 1j*rng.standard_normal(sy.nvar))*0.7
A, B, H = sy.unpack(x)
M = len(S); n = W + M
# realise u_i, v_i in C^n with E^T u = a_i, E^T v = b_i, v_i^T u_j = H_ij/2 (i<j), v_i^T u_i = 1
U = np.vstack([A, np.eye(M)])            # n x M
# V^T U = B^T A + X^T  -> want G with diag 1, upper H/2
G = np.tril(rng.standard_normal((M, M)), -1) + np.eye(M) + np.triu(H/2, 1)
Xt = G - B.T @ A                          # M x M  = X^T
V = np.vstack([B, Xt.T])
assert np.allclose(V.T @ U, G)
vals, tr = quad_forms(m)
err = 0
for f in range(vals.shape[0]):
    eps = vals[f, S]
    P = np.eye(n, dtype=complex)
    for i in range(M):                    # ordered product, i = 0 leftmost
        if eps[i]:
            P = P @ (np.eye(n) - 2*np.outer(U[:, i], V[:, i]))
    comp = P[:W, :W]
    D = np.nonzero(eps)[0]
    if len(D):
        F = A[:, D] @ np.linalg.inv(np.eye(len(D)) + H[np.ix_(D, D)]) @ B[:, D].T
    else:
        F = np.zeros((W, W))
    err = max(err, np.abs(comp - (np.eye(W) - 2*F)).max())
print("formula vs product, max err:", err)
# Jacobian check
r, J = sy.residual(x, jac=True)
d = (rng.standard_normal(sy.nvar) + 1j*rng.standard_normal(sy.nvar))
h = 1e-6
r2 = sy.residual(x + h*d)
print("jacobian check:", np.abs((r2 - r)/h - J @ d).max(), np.abs(J@d).max())
# plain program: W=2, S = [1,2,4,1,2,4]: array 0 uses moves 0..2, array 1 uses 3..5
S2 = [1, 2, 4, 1, 2, 4]; sy2 = System(3, S2, 2)
A = np.zeros((2, 6), complex); B = np.zeros((2, 6), complex); H = np.zeros((6, 6), complex)
for l, idx in enumerate([[0, 1, 2], [3, 4, 5]]):
    for i in idx: A[l, i] = 1; B[l, i] = 1
    for a in range(3):
        for b in range(a+1, 3): H[idx[a], idx[b]] = 2   # v_i^T u_j = 1 on one register -> H = 2G = 2
x2 = np.concatenate([A.ravel(), B.ravel(), H[np.triu_indices(6, 1)]])
print("plain residual:", np.abs(sy2.residual(x2)).max())
print("deletion bounds: plain2 block3", deletion_bound(3, S2), " U3", deletion_bound(3, S), " U4", deletion_bound(4, list(range(1, 16))), " U5", deletion_bound(5, list(range(1, 32))))
