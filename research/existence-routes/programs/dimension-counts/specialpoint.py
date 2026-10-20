# Is a saving solution in special position?  Block 2, rep VV (in-rep Theorem R allows 1.5W, so 2 in 3 is solvable
# there; meta-border found it exactly).  Find an exact 2-in-3 solution with meta-border's optimiser (bsearch.py),
# convert it to dimcount's coordinates, and compare the Jacobian rank AT the solution with the generic rank.
import sys, os, json, numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.join(HERE, '..', 'border'))
import bsearch
from dimcount import Shape, rank_sv
W, T, rep = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else 'VV'
bsearch.algebra_basis.__defaults__  # (bsearch caches algebras to files in the cwd; run from a scratch dir)
Pb = bsearch.Problem(2, rep, W, T)
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
best = None
for start in range(6):
    x = rng.normal(size=Pb.npar) * 0.3
    res = minimize(Pb.loss_grad, x, jac=True, method='L-BFGS-B', options={'maxiter': 6000, 'maxfun': 24000, 'gtol': 1e-16, 'ftol': 1e-30})
    if best is None or res.fun < best.fun: best = res
    if res.fun < 1e-24: break
print('bsearch loss', best.fun, 'max coef', float(np.abs(best.x).max()))
Mb = Pb.mats(Pb.coeffs(best.x))          # slot matrices, same slot order as dimcount
S = Shape(2, rep, 'engine', W, [1] * T, [2] * W, None, 0)
coef = []; k = 0
for nsl in S.slots:
    coef.append(np.array([np.einsum('aij,ij->a', S.A.conj(), Mb[k + i]) for i in range(nsl)])); k += nsl
gens, outs = S.forward(coef)
res = max(float(np.abs(outs[t] - S.target(t)).max()) for t in range(W))
print('residual in dimcount coordinates', res)
Eeff = W * W * S.gU
J0 = S.jac_sketch(coef, Eeff + 24); r0, g0 = rank_sv(J0)
Jg = S.jac_sketch(S.random_point(), Eeff + 24); rg, gg = rank_sv(Jg)
print(json.dumps({'rep': rep, 'shape': f'{W} in {T}', 'loss': float(best.fun), 'residual': res, 'Eeff': Eeff,
                  'rank_at_solution': r0, 'gap_solution': g0, 'generic_rank': rg, 'gap_generic': gg}))
