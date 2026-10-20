# Solve the commutative-slice transfer system at ONE linear part lambda (e3-slice6).
# With free shifts and copies, every C[Sh]-linear map is free, so the slice problem splits into 2^m independent
# systems, one per lambda in F_2^m.  At lambda the characters are psi = [lambda.z] + 2Q(z) (Q binary quadratic):
#   w_i = i^[lambda.z_i] (-1)^Q(z_i),   I + A C (I - (Gam/2) C)^-1 B = i^|lambda| (-1)^trQ I,   C = diag(w - 1).
# Usage: lamsolve.py m W dirs(comma) lambda starts seed
import sys
import numpy as np
from scipy.optimize import least_squares
m, W = int(sys.argv[1]), int(sys.argv[2])
dirs = [int(t) for t in sys.argv[3].split(",")]
lam = int(sys.argv[4]); starts = int(sys.argv[5]); seed = int(sys.argv[6])
M = len(dirs)
dot = lambda a, b: bin(a & b).count("1") & 1
pairs = [(j, k) for j in range(m) for k in range(j, m)]
rows = []
for mask in range(1 << len(pairs)):
    q = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
    Qz = [sum((z >> j & 1) & (z >> k & 1) for (j, k) in q) % 2 for z in dirs]
    trq = sum(1 for (j, k) in q if j == k) % 2
    w = np.array([(1j ** dot(lam, z)) * (-1) ** Qz[i] for i, z in enumerate(dirs)])
    target = (1j ** bin(lam).count("1")) * (-1) ** trq
    rows.append((w - 1, target))
low = [(j, i) for i in range(M) for j in range(i + 1, M)]
nc = W * M + M * W + len(low)
def unpack(v):
    v = v[:nc] + 1j * v[nc:]
    A = v[:W * M].reshape(W, M); B = v[W * M:2 * W * M].reshape(M, W)
    G = np.zeros((M, M), complex)
    for t, (j, i) in enumerate(low): G[j, i] = v[2 * W * M + t]
    return A, B, G
def res(v):
    A, B, G = unpack(v)
    out = []
    for c, tg in rows:
        C = np.diag(c)
        P = np.eye(W) + A @ C @ np.linalg.solve(np.eye(M) - 0.5 * G @ C, B)
        out.append((P - tg * np.eye(W)).ravel())
    r = np.concatenate(out); return np.concatenate([r.real, r.imag])
rng = np.random.default_rng(seed)
best = None
for st in range(starts):
    v0 = rng.standard_normal(2 * nc) * 0.7
    s = least_squares(res, v0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    c = np.linalg.norm(s.fun)
    if best is None or c < best[0]: best = (c, s.x)
    if c < 1e-10: break
print(f"m={m} W={W} dirs={dirs} lambda={lam:0{m}b} best residual {best[0]:.3e} after {st+1} starts", flush=True)
np.save(f"lam_sol_m{m}_W{W}_l{lam}_{'-'.join(map(str,dirs))}.npy", best[1])
