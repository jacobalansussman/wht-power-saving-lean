#!/usr/bin/python3
"""joint/cert/staticlp.py -- agent "cert".  The static inequalities  sum_w n_w c_w(U) >= W D(U)  (one per set U of
blocks inside one family) as a system in the numbers n_w of moves per direction.
(1) Exact rational linear programme: least sum of n_w over real n_w >= 0 at W = 1 (the bound then scales with W),
    with the dual solution = an explicit combination of inequalities (a hand proof of k >= rho * W).
(2) Whole numbers: for W = 1..12 the least K for which some multiset of K directions satisfies every static inequality."""
import sys, itertools
from fractions import Fraction as Fr
v = [list(map(int, l.split())) for l in open("model_cert.txt")]
cap = {w: v[15 * (w - 1):15 * w] for w in range(1, 8)}; dem = v[105:120]
def nm(b): return "(0,0)" if b == 0 else ("(%d,0)" % b if b < 8 else "(%d,1)" % (b - 7))
ineq = {}
for comp in (list(range(8)), list(range(8, 15))):
    for r in range(1, len(comp)):
        for U in itertools.combinations(comp, r):
            out = [s for s in comp if s not in U]
            D = sum(dem[s][d] for s in out for d in U); c = tuple(sum(cap[w][s][d] for s in out for d in U) for w in range(1, 8))
            ineq.setdefault((c, D), []).append(U)
I = sorted(ineq)
print("distinct static inequalities: %d" % len(I))
# (2) whole numbers
def ok(n, W): return all(sum(n[i] * c[i] for i in range(7)) >= W * D for c, D in I)
def comps(K, parts):
    if parts == 1: yield (K,); return
    for a in range(K + 1):
        for rest in comps(K - a, parts - 1): yield (a,) + rest
maxW = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for W in range(1, maxW + 1):
    K = 0
    while True:
        sols = [n for n in comps(K, 7) if ok(n, W)]
        if sols: break
        K += 1
    print("W = %2d: least K with a multiset passing every static cut: %2d  (plain schedule: %2d ; a saving of 1e-1 needs K <= %2d) ; multisets passing at that K: %d, e.g. n_1..n_7 = %s" % (W, K, 3 * W, int(2.6883 * W), len(sols), sols[0]), flush=True)
# (1) linear programme by scipy, then the dual made rational and verified exactly
try:
    import numpy as np
    from scipy.optimize import linprog
    A = np.array([[-x for x in c] for c, D in I], float); b = np.array([-D for c, D in I], float)
    r = linprog(np.ones(7), A_ub=A, b_ub=b, bounds=[(0, None)] * 7, method="highs")
    print("LP at W = 1: least sum of n_w = %.6f ; n = %s" % (r.fun, [round(float(x), 4) for x in r.x]))
    lam = [-float(m) for m in r.ineqlin.marginals]
    lamq = [Fr(l).limit_denominator(1000) for l in lam]
    lhs = [sum(lamq[i] * I[i][0][w] for i in range(len(I))) for w in range(7)]
    rhs = sum(lamq[i] * I[i][1] for i in range(len(I)))
    print("dual, made rational: coefficients of n_1..n_7 after combining = %s (each must be <= 1) ; right side = %s * W" % ([str(x) for x in lhs], rhs))
    print("exact check of the combination: %s  ->  K = sum n_w >= %s * W = %.6f * W" % (all(x <= 1 for x in lhs) and all(l >= 0 for l in lamq), rhs, float(rhs)))
    for i in range(len(I)):
        if lamq[i] != 0: print("   weight %s on: sum_w n_w * %s >= %d * W   (U = %s)" % (lamq[i], list(I[i][0]), I[i][1], " ".join(nm(x) for x in ineq[I[i]][0])))
except Exception as e:
    print("LP part failed:", repr(e))
