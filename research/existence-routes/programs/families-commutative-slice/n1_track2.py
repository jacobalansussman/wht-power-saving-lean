# N1 at block m: follow LM restarts and report a GAUGE-INVARIANT size: max over all W-subsets I of
# |det A_I det B_I| (invariant under A -> gA, B -> hB with det g det h = 1), to tell a border from gauge drift.
import sys, itertools, numpy as np, warnings
warnings.filterwarnings("ignore")
from n1 import build, resid
m, r, seed, nst, iters = map(int, sys.argv[1:6])
M, W, idx, phi = build(m, r)
subs = np.array(list(itertools.combinations(range(M), W)))
def inv_size(x):
    A = x[:W*M].reshape(W, M); B = x[W*M:].reshape(W, M)
    da = np.linalg.det(np.transpose(A[:, subs], (1, 0, 2))); db = np.linalg.det(np.transpose(B[:, subs], (1, 0, 2)))
    return np.abs(da*db).max()
rng = np.random.default_rng(seed)
for k in range(nst):
    x = (rng.standard_normal(2*W*M) + 1j*rng.standard_normal(2*W*M))/np.sqrt(2)
    mu = 1e-3; rr, J = resid(x, M, W, idx, phi, True); f = np.vdot(rr, rr).real
    log = []
    for it in range(iters):
        JH = J.conj().T; JJ = JH @ J; g = JH @ rr; dg = np.real(np.diag(JJ))
        ok = False
        while mu < 1e12:
            xn = x + np.linalg.solve(JJ + mu*np.diag(dg + 1e-12), -g)
            rn = resid(xn, M, W, idx, phi); fn = np.vdot(rn, rn).real
            if fn < f: x, f, ok = xn, fn, True; mu = max(mu/3, 1e-12); break
            mu *= 4
        if not ok or f < 1e-26: break
        rr, J = resid(x, M, W, idx, phi, True)
        if it % (iters//6) == 0: log.append((it, float('%.1e' % f), float('%.2e' % inv_size(x))))
    print(f"m={m} r={r} seed {seed} start {k}: |res|^2 {f:.3e} after {it} its; inv size {inv_size(x):.3e}; track {log}", flush=True)
    if f < 1e-20: np.save(f"n1_sol_m{m}_r{r}_s{seed}_{k}.npy", x)
