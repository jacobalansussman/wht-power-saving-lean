# check Theorem 10's three facts on solutions the solver actually finds (refl.py convention: F_Q = A_D (I+H_D)^{-1} B_D^T,
# H strictly upper, K = B^T A, E = I + H - K):
#  (a) E = 0 at positions joined by a tight form; (b) E_D nonsingular for even Q; (c) rank E_D = |D| - W for odd Q.
import numpy as np, warnings
warnings.filterwarnings("ignore")
from refl2 import System2, random_x
from refl import quad_forms
cases = [(3, [1, 2, 3, 4, 5, 6, 7], 2), (3, [1, 2, 4, 3, 6, 7, 5], 2), (4, [1, 7, 2, 11, 4, 13, 8, 14], 2),
         (3, [1, 1, 2, 2, 4, 4], 2), (4, [7, 11, 13, 14], 1), (3, [1, 2, 4], 1)]
for (m, S, W) in cases:
    sy = System2(m, S, W); rng = np.random.default_rng(5)
    for k in range(8):
        x, f, it = sy.solve(random_x(sy, rng), iters=600)
        if f < 1e-20: break
    A, B, H = sy.unpack(x); K = B.T @ A; M = len(S)
    E = np.eye(M) + H - K
    vals, tr = quad_forms(m); sub = vals[:, S]; wt = sub.sum(1)
    joined = np.zeros((M, M), bool)
    for q in np.nonzero((tr == 1) & (wt == W))[0]:
        D = np.nonzero(sub[q])[0]; joined[np.ix_(D, D)] = True
    a = np.abs(E[joined]).max() if joined.any() else 0.0
    smin_even = np.inf; rank_err = 0.0; n_even = n_odd = 0
    for q in range(vals.shape[0]):
        D = np.nonzero(sub[q])[0]
        if len(D) == 0: continue
        sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
        if tr[q] == 0:
            n_even += 1; smin_even = min(smin_even, sv[-1])
        else:
            n_odd += 1; r = len(D) - W
            # singular values beyond r must vanish, the r-th must not
            tail = sv[r:].max() if r < len(sv) else 0.0
            head = sv[r - 1] if r >= 1 else np.inf
            rank_err = max(rank_err, tail / max(1.0, sv[0]))
            if head < 1e-6: print("  odd form with rank below |D|-W:", [S[i] for i in D], sv)
    print(f"m={m} S={S} W={W}: residual {f:.1e}; (a) max |E| on joined positions {a:.1e}; "
          f"(b) min singular value of E_D over {n_even} even forms {smin_even:.2e}; "
          f"(c) max relative singular value beyond |D|-W over {n_odd} odd forms {rank_err:.1e}", flush=True)
