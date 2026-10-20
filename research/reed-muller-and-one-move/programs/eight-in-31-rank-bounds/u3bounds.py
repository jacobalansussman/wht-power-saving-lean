"""u3bounds.py (r3-classes): soundness check of starrec6's bounds on GENUINE programs.  Numerical U_3 2-in-7 units (block 3,
targets I, I, all 7 directions once; ex2-families solver) give a true E.  Every bound starrec6 prints (DUMPB=1, file
fb_u3.txt: 'FB form setmask u ell') must satisfy ell <= rank E_true[S,S] <= u."""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '../shared')
from refl2 import System2, random_x
FB = sorted(set((int(t[2]), int(t[3]), int(t[4])) for t in (l.split() for l in open(sys.argv[1])) if t and t[0] == 'FB'))
m, S, W = 3, [1, 2, 3, 4, 5, 6, 7], 2
for seed in (5, 6, 7):
    sy = System2(m, S, W); rng = np.random.default_rng(seed)
    for k in range(8):
        x, f, it = sy.solve(random_x(sy, rng), iters=600)
        if f < 1e-20: break
    A, B, H = sy.unpack(x); K = B.T @ A; E = np.eye(7) + H - K
    bad = 0; tight_u = tight_l = 0
    for (mask, u, l) in FB:
        D = [S.index(z) for z in range(1, 8) if mask >> (z - 1) & 1]
        sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
        rk = int((sv > 1e-7 * max(1, sv[0])).sum())
        if not (l <= rk <= u): bad += 1; print("VIOLATION set", D, "u", u, "ell", l, "true rank", rk)
        tight_u += rk == u; tight_l += (rk == l and l > 0)
    print("seed", seed, "residual %.1e" % f, "| distinct bounds", len(FB), "violations", bad, "| u tight", tight_u, "| ell tight (ell > 0)", tight_l)
