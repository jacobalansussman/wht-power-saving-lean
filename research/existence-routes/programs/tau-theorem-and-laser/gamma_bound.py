# gamma_bound.py (meta/ex-laser): the induction constant of PROOF.md Theorem 1.
# With cutoff L = n, the frontier schedule gives T(n) <= lv(n) + sum_s w(s) T(n-s), where w is the normalised source
# distribution and lv(n) = sum_s w(s) s^tau_c (Phi conservation, checked exactly in frontier_flow.py).
# If gamma(n) = lv(n) / (n^tau - sum_s w(s) (n-s)^tau) <= C for all n > n0 and T(n) <= C n^tau for n <= n0,
# then T(n) <= C n^tau for all n.  We print sup gamma over dyadic-ish ranges up to 10^7.
import sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from frontier_flow import flow, tau_c
for (W3, W2, M) in ((1, 4, 10), (1, 1, 4), (7, 38, 88), (1, 6, 14)):
    tc = tau_c(W3, W2, M)
    out = []
    lo = 3
    for hi in (10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
        g = 0.0; arg = None
        step = max(1, (hi - lo) // 4000)
        for n in range(lo, hi + 1, step):
            U, src, lv, bad = flow(n, W3, W2, M)
            tot = sum(src.values())
            num = lv / tot
            den = n ** tc - sum(w / tot * (n - s) ** tc for s, w in src.items())
            if num / den > g: g, arg = num / den, n
        out.append(f"({lo},{hi}]: sup gamma {g:.4f} at n={arg}")
        lo = hi
    print(f"{W3}+{W2} in {M} (tau_c {tc:.4f}): " + "; ".join(out), flush=True)
