# perturb56.py (ex2-certnum): is the a2-gct V_56 5-in-14 relaxed scheme an OPEN or a CLOSED condition in the moves?
# Perturb the move coefficients x by a relative eps, rebuild the paid-move outputs, re-solve the (linear) output
# coefficients, and report the residual.  Open (generic) => residual stays at float level; closed => residual ~ eps.
import numpy as np, relax2, sys
exec(open('verify56.py').read().split("worst = max")[0])      # builds g, x, gens, solve
rng = np.random.default_rng(5)
def resid_at(xx):
    beta, gamma, Uh, nrm = g.forward(xx)
    gens2 = [g.E[l] for l in range(W)] + [Uh[i] for i in range(M)]
    return max(solve(gens2, k @ g.E[l]) for l in range(W))
print('dtype x', x.dtype, 'len', x.size, 'n', n, 'N', g.N, 'types', [(d, len(s)) for d, s in tl])
print('residual at the stored point: %.3e' % resid_at(x))
for eps in [1e-10, 1e-8, 1e-6, 1e-4, 1e-2]:
    dx = rng.normal(size=x.shape) * eps * np.abs(x).max()
    print('eps %.0e: residual %.3e' % (eps, resid_at(x + dx)), flush=True)
for d, slots in tl:
    s = slots[0]; A = np.vstack([G[s, :] for G in gens])
    sv = np.linalg.svd(A, compute_uv=False)
    print('type d=%d: generator-row matrix %s, rank %d (sv %.2e | %.2e), target rows %d' % (d, A.shape, int((sv > 1e-9 * sv[0]).sum()), sv[min(A.shape) - 1], sv[0], d))
