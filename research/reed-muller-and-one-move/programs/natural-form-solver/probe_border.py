"""probe_border.py m TSPEC ORDER seed iters: long LM run; print cost and parameter norm along the way."""
import sys, warnings, numpy as np
warnings.filterwarnings("ignore")
from inst import canonical, slice_I
from glm import Model, lm, random_start
m = int(sys.argv[1]); ks = [int(v) for v in sys.argv[2].split(",")]; order = sys.argv[3]; seed = int(sys.argv[4]); iters = int(sys.argv[5])
inst, fd = canonical(m, ks, dirs="all", order=order)
md = Model(inst.t, inst.p, inst.W, inst.M)
rng = np.random.default_rng(seed)
x = random_start(md, rng, 1, 0.3)
for blk in range(iters // 100):
    x, c, h = lm(md, x, iters=100, verbose=0)
    A, B, H = md.unpack(x)
    print(f"its {100*(blk+1)} cost {c:.3e} |A| {np.abs(A).max():.2e} |B| {np.abs(B).max():.2e} |H| {np.abs(H).max():.2e}", flush=True)
    if c < 1e-24: break
np.save(f"probe_{m}_{sys.argv[2]}_{order}_{seed}.npy", x)
