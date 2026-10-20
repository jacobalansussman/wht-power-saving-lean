# ex2-duality colab_lp/lp_big.py: ex-laser's fluid LP (grouped_lp.py copied unchanged) at larger n, for
# (i) the recycled rank-1 rule {3:1,2:2} <= {1:6} (= ex-laser's grouped test rule, LP-identical; Theorem 1 => bounded ratio)
# (ii) the genuinely grouped mixed rule {4:1,3:2} <= {1:5,2:2} (duality prediction: exponent tau_c, ratio n^o(1)).
import sys, time
sys.path.insert(0, ".")
from math import log
from grouped_lp import solve, tau_c
Ns = [int(x) for x in sys.argv[1:]] or [480, 960, 1920]
for name, Wt, P in [("recycled rank-1 {3:1,2:2} <= {1:6}", {3: 1, 2: 2}, {1: 6}),
                    ("grouped mixed {4:1,3:2} <= {1:5,2:2}", {4: 1, 3: 2}, {1: 5, 2: 2})]:
    tc = tau_c(Wt, P)
    out = []
    for N in Ns:
        t0 = time.time(); c = solve(N, Wt, P); out.append((N, c))
        print(f"{name}: n={N} cost {c:.6f} cost/n^tau_c {c / N ** tc:.4f} ({time.time() - t0:.0f} s)", flush=True)
    print(f"{name}: tau_c {tc:.4f}; local slopes "
          f"{[round(log(out[i+1][1] / out[i][1]) / log(out[i+1][0] / out[i][0]), 4) for i in range(len(out) - 1)]}", flush=True)
