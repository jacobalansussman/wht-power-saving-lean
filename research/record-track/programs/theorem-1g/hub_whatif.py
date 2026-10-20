# hub_whatif.py (rl-prove): mixed-target price of the hub's T74 (#1839) profile and #1042 IF every idle child is a
# head/tail block on its own array (Lemma 6 hypothesis, NOT verified here) -- [what-if]; plus the constants of
# Theorem 1G (theta, eta, kappa) for these rules.  Control: moment root reproduces #1839 b = 1.6797549561e-3.
from decimal import Decimal as D, getcontext
from math import log
getcontext().prec = 60
def root(f, lo=D('0.9'), hi=D(1)):
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2
def pw(x, t): return (D(x).ln() * t).exp()
cases = {
 'T74 #1839': (73, 66304, {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328}, {48: 22748}),
 '#1042': (73, 14336, {r: 144 * n for r, n in {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}.items()},
           {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}),
}
for name, (m, W, calls, idle) in cases.items():
    f_old = lambda t: sum(n * pw(r, t) for r, n in calls.items()) + sum(n * pw(r, t) for r, n in idle.items()) - W * pw(m, t)
    ni = sum(idle.values())
    f_new = lambda t: sum(n * pw(r, t) for r, n in calls.items()) - (sum(n * pw(m - r, t) for r, n in idle.items()) + (W - ni) * pw(m, t))
    t_old = root(lambda x: -f_old(x)); t_new = root(lambda x: -f_new(x))
    kmin = min([m - r for r in idle] + [m]); rmax = max(calls)
    tau = float(t_new); theta = log(rmax) / log(kmin)
    eta = min((1 - (1 - t) ** tau) / t ** tau for t in [1 / m + (1 - 1 / m) * i / 20000 for i in range(20001)])
    kappa = -(-log(2 / eta) // log(1 / theta))
    print(f"{name}: moment root (hub engine form) saving {1 - t_old:.10e}; mixed-target what-if saving {1 - t_new:.10e} "
          f"(+{(t_old - t_new) / (1 - t_old) * 100:.2f} %); idle children {ni} on {W} arrays; separation r_max={rmax} < k_min={kmin}; "
          f"theta={theta:.4f}, eta={eta:.4f}, kappa={int(kappa)}")
