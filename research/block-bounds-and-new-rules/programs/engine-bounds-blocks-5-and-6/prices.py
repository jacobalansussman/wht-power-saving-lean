# pf-induction: engine-price caps (w-engine.md 4.1: a block of rank j costs (j/m)^z) for the block-4 hook carried
# up by projection; cap on the saving = 1 - min_j ln(f_m/f_j)/ln(m/j)  (as in b56-compute/pricecap.py).
# f_j of the projected family are certified lower bounds for j >= 2 and exact for j = 1; the price cap needs f_j
# from ABOVE for j >= 2, so these caps are only indicative (exact if the mod-p ranks are exact). Usage: python3 -I prices.py
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nrs_closed import N_closed
def step(a, m):
    return [0] + [sum(c * a[s] for s, c in N_closed(m, r).items() if 1 <= s <= m - 1) for r in range(1, m + 1)]
a = [0, 27456, 48576, 68320, 86592]
for m in range(5, 11):
    a = step(a, m)
    z = [math.log(a[m] / a[j]) / math.log(m / j) for j in range(1, m)]
    jb = min(range(1, m), key=lambda j: z[j - 1])
    print(f"m={m}: moves-only cap {1 - z[0]:.4f}; engine-priced cap {1 - z[jb - 1]:.4f} (binding rank {jb}); "
          f"f_(m-1)/f_m = {a[m-1]/a[m]:.4f}")
