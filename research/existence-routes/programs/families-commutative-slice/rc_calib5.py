# calibration of the block-5 rank-cover LP rows (rc_lp5.py) on EXISTING block-5 units: with c fixed to an existing unit's
# counts / W, the LP (piece variables free) must be FEASIBLE.  Infeasible on an existing unit = the rows are wrong.
import sys, runpy, numpy as np, warnings
warnings.filterwarnings("ignore")
from scipy.optimize import linprog
sys.argv = ["rc_lp5.py", "8"]
import io, contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    g = runpy.run_path("rc_lp5.py")
A, b, n, nv = g["A"], np.array(g["b"]), g["n"], g["nv"]
m = 5
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
def cls(z): return tuple(((z >> j) & 1) * ((z >> k) & 1) for j in range(m) for k in range(j, m))
I = tuple(1 if j == k else 0 for j in range(m) for k in range(j, m))
pool = []
while len(pool) < 300:
    k = int(rng.integers(m, m + 5)); zs = [int(z) for z in rng.integers(1, 32, size=k)]
    s = np.zeros(len(I), dtype=int)
    for z in zs: s = (s + np.array(cls(z))) % 2
    if tuple(s.tolist()) == I: pool.append(zs)
cases = [([1, 2, 4, 8, 16], 1), ([1, 2, 4, 8, 16] * 3, 3)]
for t in range(120):
    W = int(rng.integers(1, 5))
    S = sum((pool[i] for i in rng.integers(0, len(pool), size=W)), [])
    S += [int(z) for z in rng.integers(1, 32, size=int(rng.integers(0, 4)))]
    cases.append((S, W))
bad = 0
for S, W in cases:
    cnt = np.bincount(S, minlength=32)[1:] / W
    bounds = [(float(cnt[z]), float(cnt[z])) for z in range(n)] + [(None, None)] * (nv - n)
    res = linprog(np.zeros(nv), A_ub=A, b_ub=b, bounds=bounds, method="highs")
    if res.status != 0: bad += 1; print("INFEASIBLE on an existing unit:", W, sorted(S), res.status, flush=True)
print(f"{len(cases)} existing block-5 units: rank-cover LP infeasible on {bad}", flush=True)
for r, W in [(1, 7), (1, 8), (5, 36), (5, 37)]:
    cnt = np.full(n, r / W)
    bounds = [(float(v), float(v)) for v in cnt] + [(None, None)] * (nv - n)
    res = linprog(np.zeros(nv), A_ub=A, b_ub=b, bounds=bounds, method="highs")
    print(f"U_5^({r}) with W={W}: {'feasible' if res.status == 0 else 'INFEASIBLE (excluded)'}")
