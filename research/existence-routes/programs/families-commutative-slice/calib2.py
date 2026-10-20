# calibration part 2: a8 (W = 4 at plain cost, block 4) on random subsets of forms, and the infeasible controls.
import sys, numpy as np, time, warnings
warnings.filterwarnings("ignore")
from refl2 import System2, random_x
def run(m, S, W, restarts, seed=1, iters=400, label="", subset=None, tmax=None):
    rng = np.random.default_rng(seed); sy = System2(m, S, W, subset=subset, rng=np.random.default_rng(seed+100)); fs = []; t0 = time.time()
    for k in range(restarts):
        x, f, it = sy.solve(random_x(sy, rng), iters=iters, tmax=tmax); fs.append(f)
        if f < 1e-13 and subset is not None:
            full = System2(m, S, W); rf = full.residual(x); fs[-1] = (f, float(np.vdot(rf, rf).real))
    print(f"{label} m={m} W={W} M={len(S)} forms={len(sy.Ds)} vars={sy.nvar} results={fs} t={time.time()-t0:.1f}s", flush=True)
which = sys.argv[1]
if which == "a8":
    run(4, [1, 7, 2, 11, 4, 13, 8, 14, 1, 7, 2, 11, 4, 13, 8, 14], 4, 3, label="a8 4in16 b4 2x2 reps interleaved (subset 250)", subset=250, tmax=200)
if which == "bad":
    for lab, m, S, W in [("b1 2in3 b2", 2, [1, 2, 3], 2), ("b2 2in5 b3", 3, [1, 2, 4, 3, 7], 2),
                         ("b2' 2in5 b3", 3, [1, 2, 3, 4, 7], 2), ("b3 2in7 b4", 4, [1, 2, 4, 8, 7, 11, 13], 2),
                         ("b4 2in6 b3 nonplain", 3, [1, 2, 4, 3, 5, 6], 2), ("b5 3in8 b3", 3, [1, 2, 3, 4, 5, 6, 7, 7], 3)]:
        try: run(m, S, W, 6, label=lab, iters=600)
        except ValueError as e: print(lab, "->", e, flush=True)
