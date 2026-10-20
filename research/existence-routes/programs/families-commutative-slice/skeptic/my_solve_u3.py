# Own solver: find U_3 = 2 arrays in the 7 nonzero directions of F_2^3 (commutative slice) from (*) on all binary forms.
import numpy as np, sys
from scipy.optimize import least_squares
m, W = 3, 2
dirs = [1, 2, 3, 4, 5, 6, 7] if len(sys.argv) < 2 else [int(t) for t in sys.argv[1].split(",")]
M = len(dirs)
pairs = [(j, k) for j in range(m) for k in range(j, m)]
forms = []
for mask in range(1 << len(pairs)):
    q = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
    D = [i for i, z in enumerate(dirs) if sum((z >> j & 1) & (z >> k & 1) for (j, k) in q) % 2]
    tr = sum(1 for (j, k) in q if j == k) % 2
    forms.append((D, tr))
low = [(j, i) for i in range(M) for j in range(i + 1, M)]
nvar = 2 * (W * M + M * W + len(low))
def unpack(v):
    v = v[: nvar // 2] + 1j * v[nvar // 2:]
    A = v[: W * M].reshape(W, M); B = v[W * M: 2 * W * M].reshape(M, W)
    G = np.zeros((M, M), complex)
    for t, (j, i) in enumerate(low): G[j, i] = v[2 * W * M + t]
    return A, B, G
def res(v):
    A, B, G = unpack(v)
    out = []
    for D, tr in forms:
        if not D:
            continue
        X = A[:, D] @ np.linalg.solve(np.eye(len(D)) + G[np.ix_(D, D)], B[D, :])
        R = X - tr * np.eye(W)
        out.append(R.ravel())
    r = np.concatenate(out)
    return np.concatenate([r.real, r.imag])
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
best = None
for start in range(40):
    v0 = rng.standard_normal(nvar) * 0.7
    try:
        s = least_squares(res, v0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    except Exception as e:
        continue
    c = np.linalg.norm(s.fun)
    if best is None or c < best[0]:
        best = (c, s.x)
    if c < 1e-11:
        break
print("dirs", dirs, "W", W, "best residual %.3e after %d starts" % (best[0], start + 1))
A, B, G = unpack(best[1])
np.savez("u3_solution.npz", A=A, B=B, Gam=G, dirs=np.array(dirs), W=W, m=m)
