# e3-boost lp.py: how deep can a relation that the known certificates allow go, with and without oversized blocks?
# For a register of m rows (sizes 1..m), net vector v = p - w (p calls >= 0, w tasks >= 0, sum w = 1), constraints
# c . v >= 0 for every applicable certificate row (certlib.certify's rows, register = m), and objective
# L(theta) = min D_theta(v).  Three classes:
#   sep   : every call rank <= every task rank (ordinary units, mixed targets, grouped blocks below the targets);
#   boost : any ranks, D(1) >= 0 (saves nothing at plain prices);
#   all   : any ranks, D(1) free.
# theta_class(m) = inf { theta : L_class(theta) < 0 } on a grid; for sep and all the negative set is an up-set.
import sys
import numpy as np
from scipy.optimize import linprog
import certlib

def rows(m, statuses):
    R, names = [], []
    unit = {d: 1 for d in range(1, m + 1)}
    for name, b, st, _ in certlib.certify(unit, m, statuses):
        if name == "frame-count":
            R.append([1.0] * m); names.append("frame-count"); continue
        P = next(P for (n, M, P, s) in certlib.PROFILES if n == name)
        R.append([float(P[b * d - 1] / P[-1]) for d in range(1, m + 1)]); names.append(f"{name}x{b}")
    R = np.array(R); R = R / np.abs(R).max(axis=1, keepdims=True)
    return R, names

def solve(m, theta, R, cls, t_sep=None):
    # variables: p_1..p_m, w_1..w_m
    phi = np.array([d ** theta for d in range(1, m + 1)])
    c = np.concatenate([phi, -phi])
    A_ub = np.hstack([-R, R])              # -(c.(p-w)) <= 0
    b_ub = np.zeros(len(R))
    if cls == "boost":
        one = np.array([float(d) for d in range(1, m + 1)])
        A_ub = np.vstack([A_ub, np.concatenate([-one, one])]); b_ub = np.append(b_ub, 0.0)
    A_eq = np.concatenate([np.zeros(m), np.ones(m)])[None, :]
    bounds = [(0, None)] * (2 * m)
    if cls == "sep":   # calls only at d <= t, tasks only at d >= t
        bounds = [(0, None) if d <= t_sep else (0, 0) for d in range(1, m + 1)] + \
                 [(0, None) if d >= t_sep else (0, 0) for d in range(1, m + 1)]
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=bounds, method="highs")
    if r.status != 0: return None, None
    p, w = r.x[:m], r.x[m:]
    return r.fun, {d + 1: round(p[d] - w[d], 6) for d in range(m) if abs(p[d] - w[d]) > 1e-9}

def L(m, theta, R, cls):
    if cls != "sep": return solve(m, theta, R, cls)
    best = (1e9, None)
    for t in range(1, m + 1):
        f, v = solve(m, theta, R, "sep", t)
        if f is not None and f < best[0]: best = (f, v)
    return best

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "rigorous"
    statuses = certlib.RIGOROUS if mode == "rigorous" else certlib.RIGOROUS + ("sampled", "conj")
    grid = [i / 1000 for i in range(500, 1001, 5)] + [0.9983, 0.999]
    for m in range(2, 9):
        R, names = rows(m, statuses)
        out = {}
        for cls in ("sep", "all", "boost"):
            neg = [t for t in grid if (L(m, t, R, cls)[0] or 0) < -1e-9]
            out[cls] = neg
        def desc(neg):
            if not neg: return "none on [0.5, 1]"
            return f"from {min(neg):.3f} ({len(neg)} grid points of {len(grid)})"
        print(f"m={m} [{mode}; {len(names)} certificate rows: {', '.join(names)}]")
        for cls in ("sep", "all", "boost"):
            print(f"   {cls:5s}: L < 0 {desc(out[cls])}")
        for t in (0.9, 0.95, 0.99, 0.9983):
            for cls in ("sep", "boost"):
                f, v = L(m, t, R, cls)
                print(f"      theta={t} {cls}: min D = {f:.5f}  v = {v}")
        sys.stdout.flush()
