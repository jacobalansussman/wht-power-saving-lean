# tstab.py (ex2-certnum): solve a single-irreducible relaxation of a2-gct (relax2.Game, mode target) from a random start,
# then test stability: nudge the target by eps and re-optimise all coefficients from the solution.
import numpy as np, relax2, sys, time
from scipy.optimize import minimize
irr, W, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
n, h, k, tl = relax2.load([irr], np.random.default_rng(12345))
g = relax2.Game(n, h, k, tl, W, M); base = [k @ g.E[l] for l in range(W)]; g.set_target(base)
rng = np.random.default_rng(seed)
best = None
for st in range(3):
    x0 = rng.standard_normal(g.npar) / np.sqrt(max(1, max(d for d, _ in tl)))
    r = minimize(g.fg, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=3000, maxcor=50, ftol=1e-30, gtol=1e-18, maxfun=9000))
    if best is None or r.fun < best.fun: best = r
    if r.fun < 1e-26: break
x = best.x
line = 'irr %d (n=%d) %d in %d: solved to %.1e (|x| %.1f);' % (irr, n, W, M, best.fun, np.linalg.norm(x))
for eps in [1e-6, 1e-4]:
    outs = []
    for T in base:
        P = rng.normal(size=T.shape); outs.append(T + eps * P * np.linalg.norm(T) / np.linalg.norm(P))
    g.set_target(outs)
    r = minimize(g.fg, x, jac=True, method='L-BFGS-B', options=dict(maxiter=1500, maxcor=50, ftol=1e-30, gtol=1e-18, maxfun=5000))
    line += ' nudge %.0e -> %.1e (|dx|/|x| %.1e);' % (eps, r.fun, np.linalg.norm(r.x - x) / np.linalg.norm(x))
    g.set_target(base)
print(line, flush=True)
