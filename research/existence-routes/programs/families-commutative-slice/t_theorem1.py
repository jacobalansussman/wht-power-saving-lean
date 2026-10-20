# check Theorem 1's key step on solutions the solver actually finds: for every tight odd form, K[D,D] = I + H_D
# (refl.py convention: F_Q = A_D (I + H_D)^{-1} B_D^T, H strictly upper, K = B^T A), and rank K <= W.
import numpy as np, warnings
warnings.filterwarnings("ignore")
from refl2 import System2, random_x
from refl import quad_forms
for (m, S, W) in [(3, [1, 2, 3, 4, 5, 6, 7], 2), (4, [1, 7, 2, 11, 4, 13, 8, 14], 2), (3, [1, 1, 2, 2, 4, 4], 2)]:
    sy = System2(m, S, W); rng = np.random.default_rng(5)
    for k in range(5):
        x, f, it = sy.solve(random_x(sy, rng), iters=600)
        if f < 1e-20: break
    A, B, H = sy.unpack(x); K = B.T @ A
    vals, tr = quad_forms(m); sub = vals[:, S]; wt = sub.sum(1)
    errs = [np.abs(K[np.ix_(D, D)] - (np.eye(len(D)) + H[np.ix_(D, D)])).max()
            for D in (np.nonzero(sub[q])[0] for q in np.nonzero((tr == 1) & (wt == W))[0])]
    sv = np.linalg.svd(K, compute_uv=False)
    print(f"m={m} S={S} W={W}: residual {f:.1e}; tight forms {len(errs)}; max |K[D,D] - (I+H_D)| = {max(errs):.1e}; "
          f"singular values of K beyond W: {sv[W:].max() if len(sv) > W else 0:.1e}")
