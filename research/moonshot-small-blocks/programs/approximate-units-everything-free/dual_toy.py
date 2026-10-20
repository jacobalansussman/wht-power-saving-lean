"""dual_toy.py -- w7-xwborder.  Proposition Lambda0: over the dual numbers C[eps]/(eps^2), the two layer equations that
give Theorem 2W (B G = 0 and B n G = eps K, n strictly lower triangular) admit a solution with M = 3, W = 2 (K = 1).
Exact integer arithmetic.  Also: the same equations with eps = 0 (exact case) have no solution with M = 3 (rank
argument), checked here by a random search of the rank condition as a sanity control."""
import numpy as np
e = np.eye(3, dtype=np.int64); f2 = np.eye(2, dtype=np.int64)
b0 = np.outer(f2[1], e[2]); b1 = np.outer(f2[0], e[1])          # 2 x 3
g0 = np.outer(e[0], f2[0]); g1 = np.outer(e[1], f2[1])          # 3 x 2
n0 = np.zeros((3, 3), np.int64); n0[1, 0] = 1; n0[2, 1] = 1; n1 = np.zeros((3, 3), np.int64)
def mul(x, y):            # (x0 + eps x1)(y0 + eps y1) mod eps^2
    return (x[0] @ y[0], x[0] @ y[1] + x[1] @ y[0])
B, G, N = (b0, b1), (g0, g1), (n0, n1)
BG = mul(B, G); BNG = mul(mul(B, N), G)
print("B G   order 0:", BG[0].tolist(), " order 1:", BG[1].tolist())
print("B n G order 0:", BNG[0].tolist(), " order 1:", BNG[1].tolist(), "(want 0 and I)")
ok = (not BG[0].any()) and (not BG[1].any()) and (not BNG[0].any()) and (BNG[1] == f2).all()
print("Proposition Lambda0 example:", "HOLDS" if ok else "FAILS")
# control: exact case b g = 0 and b n g = I with M = 3, W = 2 is impossible (r(b)+r(g) <= 3, both >= 2):
rng = np.random.default_rng(1); best = 1e9
for t in range(20000):
    g = rng.standard_normal((3, 2)); 
    # b must kill im(g): b = c * (orthogonal complement row)
    q, _ = np.linalg.qr(np.c_[g, rng.standard_normal(3)]); w = q[:, 2]
    c = rng.standard_normal(2); b = np.outer(c, w)
    n = np.tril(rng.standard_normal((3, 3)), -1)
    best = min(best, np.linalg.svd(b @ n @ g, compute_uv=False)[1])
print("control (exact case, M=3, W=2): smallest second singular value of b n g over 20000 random tries = %.1e (rank <= 1 always)" % best)
