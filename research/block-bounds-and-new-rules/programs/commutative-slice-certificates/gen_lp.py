# pf-logging, 14 Oct 2026.  GENERATOR (untrusted) of a single "mixed strategy" certificate for g-play's commutative-slice
# test, valid for EVERY number W of arrays at once.  For block size m, functionals y in F_2^{m + C(m,2)} with
# y.Phi(kappa) = 1 give quadratic functions q_y(z) = y.Phi(lift z) on directions z != 0.  An LP finds a probability
# distribution mu on such y minimising v = max_z E_mu[q_y(z)].  Writes cert_lp_m{m}.json with exact rationals
# (limit_denominator; the checker re-verifies everything exactly).   Usage: python3 (with numpy) -I gen_lp.py m OUT.json
import sys, itertools, json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
m = int(sys.argv[1]); out = sys.argv[2]
pairs = list(itertools.combinations(range(m), 2)); N = m + len(pairs)
def phi(z):
    v = z
    for b, (j, l) in enumerate(pairs):
        if (z >> j) & 1 and (z >> l) & 1: v |= 1 << (m + b)
    return v
kap = (1 << m) - 1
par = lambda x: bin(x).count('1') & 1
ys = [y for y in range(1 << N) if par(y & kap)]
A = np.array([[par(y & phi(z)) for y in ys] for z in range(1, 1 << m)], dtype=float)
ny = len(ys); c = np.zeros(ny + 1); c[-1] = 1
r = linprog(c, A_ub=np.hstack([A, -np.ones((A.shape[0], 1))]), b_ub=np.zeros(A.shape[0]),
            A_eq=np.hstack([np.ones((1, ny)), np.zeros((1, 1))]), b_eq=[1],
            bounds=[(0, None)] * ny + [(None, None)], method='highs')
mu = {ys[i]: Fraction(r.x[i]).limit_denominator(100000) for i in range(ny) if r.x[i] > 1e-9}
tot = sum(mu.values())
mu = {y: q / tot for y, q in mu.items()}                      # renormalise exactly
v = max(sum((q for y, q in mu.items() if par(y & phi(z))), Fraction(0)) for z in range(1, 1 << m))
json.dump({"m": m, "v": str(v), "mu": {str(y): str(q) for y, q in mu.items()},
           "note": "y as an integer: bits 0..m-1 = linear part, bit m+b = pair b of itertools.combinations(range(m),2)"},
          open(out, 'w'), indent=1)
print(f"m={m}: LP value {r.fun:.9f}; exact certificate v = {v} = {float(v):.9f}; support {len(mu)}; bound n < W / v = {1/v} W")
