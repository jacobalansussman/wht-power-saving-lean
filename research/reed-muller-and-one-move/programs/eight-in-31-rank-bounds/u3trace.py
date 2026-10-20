"""u3trace.py: in genuine numerical U_3 2-in-7 units (ex2-families solver), print diag(E) and per-array traces a_l.b_l."""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '../shared')
from refl2 import System2, random_x
m, S, W = 3, [1, 2, 3, 4, 5, 6, 7], 2
for seed in (5, 6, 7, 8, 9):
    sy = System2(m, S, W); rng = np.random.default_rng(seed)
    for k in range(8):
        x, f, it = sy.solve(random_x(sy, rng), iters=600)
        if f < 1e-20: break
    A, B, H = sy.unpack(x); K = B.T @ A; E = np.eye(7) + H - K
    kap = B * A
    print("seed", seed, "res %.1e" % f, "diagE", np.round(np.diag(E), 6), "traces a_l.b_l", np.round(kap.sum(axis=1), 6), "trK %.6f" % np.trace(K))
    print("   kappa rows", np.round(kap, 4).tolist())
