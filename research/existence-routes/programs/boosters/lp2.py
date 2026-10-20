# e3-boost lp2.py: (1) boosters on the whole of [0, 1]; (2) boosters whose ranks stay below the register size
# (spare rows used only by free steps); (3) the COVER test: one relation, allowed by the certificates, with D(1) >= 0
# and D < 0 on all of [theta0, 0.9983] (with the known units, which beat every theta > 0.99833, such a relation would
# push the exponent down to theta0: Corollary 6 of ex2-duality); smallest theta0 by bisection.
import sys
import numpy as np
from scipy.optimize import linprog
import certlib
from lp import rows

def boost_min(m, theta, R, rmax=None):
    rmax = rmax or m
    phi = np.array([d ** theta for d in range(1, m + 1)]); one = np.arange(1, m + 1, dtype=float)
    c = np.concatenate([phi, -phi])
    A_ub = np.vstack([np.hstack([-R, R]), np.concatenate([-one, one])]); b_ub = np.zeros(len(R) + 1)
    A_eq = np.concatenate([np.zeros(m), np.ones(m)])[None, :]
    bounds = [(0, None) if d <= rmax else (0, 0) for d in range(1, m + 1)] * 2
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=bounds, method="highs")
    return r.fun, (r.x[:m] - r.x[m:])

def cover(m, R, th0, th1=0.9983, n=60):
    """max eps s.t. D(theta_i) <= -eps on a grid of [th0, th1]; D(1) >= 0; certificates; sum w = 1 (tasks normalised)."""
    ths = [th0 + (th1 - th0) * i / n for i in range(n + 1)]
    one = np.arange(1, m + 1, dtype=float)
    A_ub, b_ub = [], []
    for row in R: A_ub.append(np.concatenate([-row, row, [0.0]])); b_ub.append(0.0)      # certificates
    A_ub.append(np.concatenate([-one, one, [0.0]])); b_ub.append(0.0)                    # D(1) >= 0
    for t in ths:
        phi = one ** t
        A_ub.append(np.concatenate([phi, -phi, [1.0]])); b_ub.append(0.0)                # D(t) + eps <= 0
    A_eq = np.concatenate([np.zeros(m), np.ones(m), [0.0]])[None, :]
    c = np.zeros(2 * m + 1); c[-1] = -1.0
    bounds = [(0, None)] * (2 * m) + [(None, 1.0)]
    r = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=[1.0], bounds=bounds, method="highs")
    if r.status != 0: return None, None
    return -r.fun, (r.x[:m] - r.x[m:2 * m])

def fmt(v):
    return {d + 1: round(float(x), 5) for d, x in enumerate(v) if abs(x) > 1e-7}

if __name__ == "__main__":
    grid = [i / 200 for i in range(1, 200)] + [0.999]
    modes = sys.argv[1:] or ["rigorous", "argued", "all"]
    for mode in modes:
        statuses = {"rigorous": certlib.RIGOROUS, "argued": certlib.RIGOROUS + ("argued",),
                    "all": certlib.RIGOROUS + ("argued", "sampled", "conj")}[mode]
        print(f"===== certificates: {mode} =====")
        for m in range(2, 9):
            R, names = rows(m, statuses)
            for rmax in range(2, m + 1):
                neg = [t for t in grid if boost_min(m, t, R, rmax)[0] < -1e-9]
                if rmax == m or neg:
                    print(f"m={m} ranks<={rmax}: boosters negative at {len(neg)}/{len(grid)} grid points of (0,1)"
                          + (f", from {min(neg):.3f} to {max(neg):.3f}" if neg else ""))
            # cover: smallest theta0 with eps > 0
            e, v = cover(m, R, 0.9)
            lo, hi = 0.0, 0.9983
            if e is not None and e > 1e-9:
                hi = 0.9
                lo = 0.0
            else:
                lo = 0.9
            for _ in range(30):
                mid = (lo + hi) / 2
                e2, _ = cover(m, R, mid)
                if e2 is not None and e2 > 1e-9: hi = mid
                else: lo = mid
            e3, v3 = cover(m, R, hi)
            print(f"   cover [0.9, 0.9983]: eps = {e:.3g} (eps > 0: one certificate-allowed booster negative on the whole "
                  f"window); smallest theta0 with a cover of [theta0, 0.9983]: {hi:.4f} (eps {e3:.3g}), v = {fmt(v3)}")
            sys.stdout.flush()
