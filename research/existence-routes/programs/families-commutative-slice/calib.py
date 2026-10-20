# calibration of the reflection solver (refl.py): known-feasible and known-infeasible commutative-slice targets.
import sys
from refl2 import System2, random_x
import numpy as np, time, warnings
warnings.filterwarnings("ignore")
def run(m, S, W, restarts=6, seed=1, iters=400, label="", subset=None):
    rng = np.random.default_rng(seed); sy = System2(m, S, W, subset=subset); fs = []; t0 = time.time()
    for k in range(restarts):
        x, f, it = sy.solve(random_x(sy, rng), iters=iters); fs.append(f)
    print(f"{label} m={m} W={W} M={len(S)} S={S} forms={len(sy.Ds)} vars={sy.nvar} best={min(fs):.3e} all={[float('%.2e' % v) for v in sorted(fs)]} t={time.time()-t0:.1f}s", flush=True)
cases_ok = [
    ("a1 plain 1in2 b2", 2, [1, 2], 1),
    ("a2 plain 1in3 b3", 3, [1, 2, 4], 1),
    ("a3 1in4 (7,11,13,14)", 4, [7, 11, 13, 14], 1),
    ("a4 plain 2in4 b2 interleaved", 2, [1, 2, 1, 2], 2),
    ("a4' plain 2in4 b2 [1,1,2,2]", 2, [1, 1, 2, 2], 2),
    ("a5 plain 2in6 b3 interleaved", 3, [1, 1, 2, 2, 4, 4], 2),
    ("a6 U3 2in7 lex", 3, [1, 2, 3, 4, 5, 6, 7], 2),
    ("a6' U3 2in7 singer", 3, [1, 2, 4, 3, 6, 7, 5], 2),
    ("a7 2in8 b4 two reps interleaved", 4, [1, 7, 2, 11, 4, 13, 8, 14], 2),
    ("a8 4in16 b4 2x2 reps interleaved", 4, [1, 7, 2, 11, 4, 13, 8, 14, 1, 7, 2, 11, 4, 13, 8, 14], 4),
]
cases_bad = [
    ("b1 2in3 b2", 2, [1, 2, 3], 2),
    ("b2 2in5 b3", 3, [1, 2, 4, 3, 7], 2),
    ("b2' 2in5 b3 (7 dirs minus 2)", 3, [1, 2, 3, 4, 7], 2),
    ("b3 2in7 b4", 4, [1, 2, 4, 8, 7, 11, 13], 2),
    ("b4 2in6 b3 nonplain", 3, [1, 2, 4, 3, 5, 6], 2),
    ("b5 3in8 b3", 3, [1, 2, 3, 4, 5, 6, 7, 7], 3),
]
which = sys.argv[1] if len(sys.argv) > 1 else "all"
R = int(sys.argv[2]) if len(sys.argv) > 2 else 6
for lab, m, S, W in (cases_ok if which in ("ok", "all") else []) + (cases_bad if which in ("bad", "all") else []):
    try:
        run(m, S, W, restarts=R, seed=1, iters=400, label=lab)
    except ValueError as e:
        print(lab, "->", e, flush=True)
