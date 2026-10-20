"""u3true.py: a genuine slice unit U_3 2 in 7 (all 7 directions once, targets I, I) found numerically by ex2-families'
solver (refl2.System2, as in t_theorem10.py); its E = I + H - K (refl convention); check (1) T10 with OWN parities
(vgen, block 3, target k = 0 = I), (2) every zero cell own vgen derives for targets (0,0) is ~0 in the true E."""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '../../shared')
from refl2 import System2, random_x
import vgen
m, S, W = 3, [1, 2, 3, 4, 5, 6, 7], 2
for seed in (5, 6, 7):
    sy = System2(m, S, W); rng = np.random.default_rng(seed)
    for k in range(8):
        x, f, it = sy.solve(random_x(sy, rng), iters=600)
        if f < 1e-20: break
    A, B, H = sy.unpack(x); K = B.T @ A; E = np.eye(7) + H - K
    bl = vgen.Block(3); w = bl.CLSPAR[:, [0, 0]].sum(axis=1); r = w[bl.VID]
    bad = 0
    for Q in range(1, bl.NF):
        D = [S.index(z) for z in range(1, 8) if int(bl.MASK[Q]) >> (z - 1) & 1]
        sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
        rk = int((sv > 1e-7 * max(1, sv[0])).sum())
        if rk != len(D) - r[Q]: bad += 1
    res = vgen.run(bl, [0, 0]); Z = vgen.LASTZ
    mx = max(abs(E[S.index(a), S.index(b)]) for (a, b) in Z)
    nz = max(abs(E[i, j]) for i in range(7) for j in range(7))
    print("seed", seed, "residual %.1e" % f, "T10 mismatches (own parities)", bad, "| vgen", res, "| max |E| on derived zero cells %.1e" % mx, "(max |E| overall %.2f)" % nz)
