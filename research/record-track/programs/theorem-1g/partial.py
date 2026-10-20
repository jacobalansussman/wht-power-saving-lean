# partial.py (rl-v-logic): [what-if] T74 mixed price when only t of the 22,748 rank-48 idle children are
# truncatable. Any remaining rank-48 paid round breaks separation (48 >= k_min = 25), so Theorem 1G needs it split
# into two rank-24 paid rounds (24 < 25). Also the log exponent: smallest kappa with theta^kappa < eta.
from decimal import Decimal as D, getcontext
from math import log
getcontext().prec = 50
pw = lambda x, t: (D(x).ln() * t).exp()
calls = {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328}; I = 22748; W73 = 43556
def root(t_idle):
    lo, hi = D('0.99'), D(1)
    for _ in range(200):
        x = (lo + hi) / 2
        lhs = sum(n * pw(r, x) for r, n in calls.items()) + (I - t_idle) * 2 * pw(24, x)
        rhs = (W73 + I - t_idle) * pw(73, x) + t_idle * pw(25, x)
        if lhs > rhs: lo = x
        else: hi = x
    return 1 - (lo + hi) / 2
for f in (0, 0.25, 0.5, 0.75, 0.9, 1.0):
    t = round(f * I); print(f"truncatable fraction {f:.2f} (t={t}): saving {root(t):.7e}  (#1839: 1.6797550e-3)")
theta = log(8) / log(25); eta = 0.9904
print(f"theta={theta:.4f}; any kappa > ln(1/eta)/ln(1/theta) = {log(1/eta)/log(1/theta):.4f} suffices in Lemma A3 (beta = eta - theta^kappa)")
