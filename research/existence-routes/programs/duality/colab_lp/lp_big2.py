# ex2-duality colab_lp/lp_big2.py: the genuinely grouped mixed rule {4:1,3:2} <= {1:5,2:2} only, at larger n.
import sys, time
sys.path.insert(0, ".")
from math import log
from grouped_lp import solve, tau_c
Wt, P = {4: 1, 3: 2}, {1: 5, 2: 2}
tc = tau_c(Wt, P); prev = None
for N in [int(x) for x in sys.argv[1:]]:
    t0 = time.time(); c = solve(N, Wt, P)
    sl = "" if prev is None else f" local slope {log(c / prev[1]) / log(N / prev[0]):.4f}"
    print(f"grouped mixed {{4:1,3:2}} <= {{1:5,2:2}} (tau_c {tc:.4f}): n={N} cost {c:.6f} cost/n^tau_c {c / N ** tc:.4f}{sl} ({time.time() - t0:.0f} s)", flush=True)
    prev = (N, c)
