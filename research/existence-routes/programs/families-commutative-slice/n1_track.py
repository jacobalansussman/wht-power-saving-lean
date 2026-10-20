# follow N1 restarts at m = 4: residual and size of the unknowns along the iterations (genuine solution or border?)
import sys, numpy as np, warnings, time
warnings.filterwarnings("ignore")
from n1 import build, resid
m, r, seed, nst = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
M, W, idx, phi = build(m, r)
rng = np.random.default_rng(seed)
for k in range(nst):
    x = (rng.standard_normal(2*W*M) + 1j*rng.standard_normal(2*W*M))/np.sqrt(2)
    mu = 1e-3; rr, J = resid(x, M, W, idx, phi, True); f = np.vdot(rr, rr).real
    log = []
    for it in range(8000):
        JH = J.conj().T; JJ = JH @ J; g = JH @ rr; dg = np.real(np.diag(JJ))
        ok = False
        while mu < 1e12:
            xn = x + np.linalg.solve(JJ + mu*np.diag(dg + 1e-12), -g)
            rn = resid(xn, M, W, idx, phi); fn = np.vdot(rn, rn).real
            if fn < f: x, f, ok = xn, fn, True; mu = max(mu/3, 1e-12); break
            mu *= 4
        if not ok or f < 1e-26: break
        rr, J = resid(x, M, W, idx, phi, True)
        if it % 1000 == 0:
            A = x[:W*M].reshape(W, M); B = x[W*M:].reshape(W, M)
            # gauge-invariant size: products |a_z| |b_z| per point and the min singular value of A
            log.append((it, f, np.abs(x).max(), np.linalg.svd(A, compute_uv=False).min()))
    A = x[:W*M].reshape(W, M); B = x[W*M:].reshape(W, M)
    print(f"seed {seed} start {k}: final |res|^2 {f:.3e} after {it} its; max|x| {np.abs(x).max():.3e}; "
          f"track {[ (a, float('%.1e'%b), float('%.1e'%c)) for (a,b,c,d) in log]}", flush=True)
    if f < 1e-20: np.save(f"n1_sol_m{m}_r{r}_s{seed}_{k}.npy", x)
