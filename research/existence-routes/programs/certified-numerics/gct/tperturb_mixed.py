# tperturb_mixed.py (ex2-certnum): openness test for a2-gct's mixed 1 + 4 in 10 in V_56 (relax2b mode mixed1):
# perturb the target by eps, re-optimise all coefficients from the stored solution, report the objective.
import numpy as np, sys, time
import relax2b as R
from scipy.optimize import minimize
eps = float(sys.argv[1]); W, M = 5, 10
n, h, k, tl = R.load([8], np.random.default_rng(12345))
g = R.Game(n, h, k, tl, W, M)
x = np.load('sol_8_W5_M10_mixed1_s1.npy')
base = [(k if l < 1 else R.K2) @ g.E[l] for l in range(W)]
g.set_target(base); f0 = g.fg(x, need_grad=False); f0 = f0[0] if isinstance(f0, tuple) else f0
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
outs = [T + eps * (lambda P: P * np.linalg.norm(T) / np.linalg.norm(P))(rng.normal(size=T.shape)) for T in base]
g.set_target(outs); f1 = g.fg(x, need_grad=False); f1 = f1[0] if isinstance(f1, tuple) else f1
t0 = time.time()
res = minimize(g.fg, x, jac=True, method='L-BFGS-B', options={'maxiter': 400, 'maxfun': 4000, 'gtol': 1e-18, 'ftol': 1e-30})
print('1+4 in 10, V_56: objective at stored solution %.2e; target eps %.0e -> %.2e at the old point -> %.2e after %d its (|dx|/|x| %.2e, %.0fs)'
      % (f0, eps, f1, res.fun, res.nit, np.linalg.norm(res.x - x) / np.linalg.norm(x), time.time() - t0), flush=True)
