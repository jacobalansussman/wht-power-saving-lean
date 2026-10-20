# reprice_hub.py (meta/ex-laser): what mixed-target (tau-theorem) pricing WOULD give the hub's #1042 profile IF each
# idle child were a per-array tail (an array's last operations, on that array alone).  [what-if; structure unverified]
# Moment root (engine today):  sum_calls n_r r^tau + sum_idle r^tau = W 73^tau.
# Mixed pricing (Theorem 1 generalised to grouped paid blocks -- NOT proved here):
#   sum_calls n_r r^tau = sum_w (73 - idle_w)^tau.
from math import log
L1042 = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
I1042 = {1: 2249, 2: 1, 3: 400, 7: 3, 16: 47, 48: 4850}
m, W = 73, 14336
calls = {r: 144 * n for r, n in L1042.items()}
def bis(f):
    lo, hi = 0.9, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return (lo + hi) / 2
# f(tau) = paid moment - target value; increasing in tau near 1 (target side has the big ranks)
f_old = lambda t: sum(n * r ** t for r, n in calls.items()) + sum(n * r ** t for r, n in I1042.items()) - W * m ** t
nidle = sum(I1042.values())
f_new = lambda t: sum(n * r ** t for r, n in calls.items()) - (sum(n * (m - r) ** t for r, n in I1042.items()) + (W - nidle) * m ** t)
for name, f in (("moment root (today)", f_old), ("mixed-target pricing (what-if)", f_new)):
    t = bis(lambda x: -f(x))
    print(f"{name}: tau = {t:.9f}, saving = {1 - t:.6e}")
print(f"idle children {nidle} on {W} arrays (assumed at most one per array)")
