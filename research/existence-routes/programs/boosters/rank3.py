# e3-boost rank3.py: the exact decomposition behind Theorem B1 (jobs of at most 3 rows cannot boost).
# For d in {1, 2, 3}:  d^t = lf(t) * 1 + lp(t) * d + lA(t) * A(d),  A = (4, 6, 7) (certificate A, any register, scaled)
#   lf = 4 - 5 2^t + 2 3^t,  lp = 1 - 3 2^t + 2 3^t,  lA = -1 + 2 2^t - 3^t.
# Checks: the identity exactly at rational sample points is not possible (irrational), so: (1) identity to 1e-14 on a
# grid; (2) each coefficient >= 0 on [0, 1] (each is a 3-term exponential sum with 2 sign changes, so at most 2 real
# roots counted with multiplicity; the known roots and end slopes then fix the sign, see NOTE section 3);
# (3) lf/lp and lA/lp are decreasing on (0, 1], which makes every Theta_rho an interval [0, tau] (no amplification
# among such relations).
from math import log
lf = lambda t: 4 - 5 * 2 ** t + 2 * 3 ** t
lp = lambda t: 1 - 3 * 2 ** t + 2 * 3 ** t
lA = lambda t: -1 + 2 * 2 ** t - 3 ** t
A = {1: 4, 2: 6, 3: 7}
err = max(abs(lf(t) + lp(t) * d + lA(t) * A[d] - d ** t) for t in [i / 1000 for i in range(1001)] for d in (1, 2, 3))
print("max identity error on the grid:", err)
print("min of lf, lp, lA on (0,1) grid:", min(lf(i / 1000) for i in range(1, 1000)), min(lp(i / 1000) for i in range(1, 1000)),
      min(lA(i / 1000) for i in range(1, 1000)))
print("end values: lf(0), lf(1) =", lf(0), lf(1), "; lp(0), lp(1) =", lp(0), lp(1), "; lA(0), lA(1) =", lA(0), lA(1))
print("slopes: lf'(1) =", -5 * 2 * log(2) + 2 * 3 * log(3), " lp'(0) =", -3 * log(2) + 2 * log(3),
      " lA'(0) =", 2 * log(2) - log(3), " lA'(1) =", 2 * 2 * log(2) - 3 * log(3))
ts = [i / 2000 for i in range(1, 2001)]
r1 = [lf(t) / lp(t) for t in ts]; r2 = [lA(t) / lp(t) for t in ts]
print("lf/lp decreasing on grid:", all(a > b for a, b in zip(r1, r1[1:])), "; lA/lp decreasing on grid:",
      all(a > b for a, b in zip(r2, r2[1:])))
