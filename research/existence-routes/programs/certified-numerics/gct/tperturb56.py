# tperturb56.py (ex2-certnum): is the V_56 5-in-14 relaxed scheme stable under a generic perturbation of the TARGET?
# Open (locally surjective program map) => after re-optimising ALL coefficients from the stored solution the residual
# returns to float level; closed (image of positive codimension, special position) => the residual stays ~ eps.
import numpy as np, relax2, sys, time
from scipy.optimize import minimize
n, h, k, tl = relax2.load([8], np.random.default_rng(12345))
W, M = 5, 14
g = relax2.Game(n, h, k, tl, W, M)
x = np.load('sol_8_W5_M14_target_s1.npy')
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
eps = float(sys.argv[1])
outs = []
for l in range(W):
    T = k @ g.E[l]; P = rng.normal(size=T.shape); outs.append(T + eps * P * np.linalg.norm(T) / np.linalg.norm(P))
g.set_target(outs)
f0 = g.fg(x, need_grad=False) if True else None
f0 = f0[0] if isinstance(f0, tuple) else f0
t0 = time.time()
res = minimize(g.fg, x, jac=True, method='L-BFGS-B', options={'maxiter': int(sys.argv[3]) if len(sys.argv) > 3 else 400,
               'maxfun': 4000, 'gtol': 1e-18, 'ftol': 1e-30})
print('eps %.0e: objective at the old solution %.3e -> after %d L-BFGS its %.3e (|dx|/|x| = %.2e, %.0fs)' %
      (eps, f0, res.nit, res.fun, np.linalg.norm(res.x - x) / np.linalg.norm(x), time.time() - t0), flush=True)
