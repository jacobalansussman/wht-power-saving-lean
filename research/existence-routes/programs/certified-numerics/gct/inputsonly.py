# inputsonly.py IRR W M [seeds] [maxiter] -- a2-gct's relaxation (relax2.Game, target mode) with every gamma forced to 0:
# all paid moves act on the inputs alone ("parallel transfers", as in the V(x)V stand-in construction).  Reports the best
# objective; 0 means an inputs-only scheme exists numerically, a floor means chaining is needed.
import sys, numpy as np, relax2, time
from scipy.optimize import minimize
irr, W, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
seeds = int(sys.argv[4]) if len(sys.argv) > 4 else 3; maxit = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
n, h, k, tl = relax2.load([irr], np.random.default_rng(12345))
g = relax2.Game(n, h, k, tl, W, M); g.set_target([k @ g.E[l] for l in range(W)])
sd = g.sdim; mask = np.zeros(g.npar, bool); o = 0
for i in range(M):
    mask[o:o + W * sd] = True; o += W * sd       # beta_i: free
    o += i * sd                                  # gamma_i,j: forced 0
assert o == g.npar
def f2(y):
    x = np.zeros(g.npar); x[mask] = y
    fv, gr = g.fg(x)
    return fv, gr[mask]
t0 = time.time()
for sd_ in range(seeds):
    rng = np.random.default_rng(sd_ + 1)
    y0 = rng.standard_normal(mask.sum()) / np.sqrt(max(d for d, _ in tl))
    r = minimize(f2, y0, jac=True, method='L-BFGS-B', options=dict(maxiter=maxit, maxcor=50, ftol=1e-30, gtol=1e-18, maxfun=3 * maxit))
    print('irr %d %d in %d, inputs only, seed %d: objective %.3e after %d its (%.0fs)' % (irr, W, M, sd_ + 1, r.fun, r.nit, time.time() - t0), flush=True)
