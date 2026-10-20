# v-e3-boost vB1.py: Theorem B1 and Corollary B2 checked with my own code (no sympy on this Mac: the identities are
# linear in 1, x = 2^t, y = 3^t, so they are checked coefficient by coefficient with Fractions, which is exact).
from fractions import Fraction as Fr
from math import log
# polynomials in (1, x, y) as dicts
lf = {"1": 4, "x": -5, "y": 2}; lp = {"1": 1, "x": -3, "y": 2}; lA = {"1": -1, "x": 2, "y": -1}
A = {1: 4, 2: 6, 3: 7}
target = {1: {"1": 1}, 2: {"x": 1}, 3: {"y": 1}}
for d in (1, 2, 3):
    s = {k: lf.get(k, 0) + d * lp.get(k, 0) + A[d] * lA.get(k, 0) for k in ("1", "x", "y")}
    print(f"B1 identity d={d}: lf + {d} lp + {A[d]} lA = {s}  (want {target[d]}):",
          all(s[k] == target[d].get(k, 0) for k in s))
print("A(d) = 8(1 - 2^-d):", [8 - Fr(8, 2 ** d) for d in (1, 2, 3)])
a, b = log(2), log(3)
f = lambda c, t: c["1"] + c["x"] * 2 ** t + c["y"] * 3 ** t
fd = lambda c, t: c["x"] * a * 2 ** t + c["y"] * b * 3 ** t
print("endpoint data: lA(0), lA(1), lA'(0) =", f(lA, 0), f(lA, 1), round(fd(lA, 0), 4),
      "| lf(0), lf(1), lf'(1) =", f(lf, 0), f(lf, 1), round(fd(lf, 1), 4),
      "| lp(0), lp'(0), lp(1) =", f(lp, 0), round(fd(lp, 0), 4), f(lp, 1))
# Wronskians: hand algebra gives W1 = -a x + b y + (a - b) x y, W2 = 7 a x - 6 b y + 4 (b - a) x y (see NOTE)
ts = [i / 1000 for i in range(1001)]
e1 = max(abs((fd(lA, t) * f(lp, t) - f(lA, t) * fd(lp, t)) - (-a * 2**t + b * 3**t + (a - b) * 6**t)) for t in ts)
e2 = max(abs((fd(lf, t) * f(lp, t) - f(lf, t) * fd(lp, t)) - (7*a * 2**t - 6*b * 3**t + 4*(b - a) * 6**t)) for t in ts)
print(f"Wronskian closed forms: max error on 1001 points {e1:.1e}, {e2:.1e}")
g = lambda t: a * (1 - 3 ** -t) - b * (1 - 2 ** -t); h = lambda t: 7*a*3**-t - 6*b*2**-t + 4*(b - a)
print("g < 0 on (0,1]:", all(g(t) < 0 for t in ts[1:]), " h < 0 on [0,1]:", all(h(t) < 0 for t in ts),
      " h(0), h(1) =", round(h(0), 4), round(h(1), 4), " h' root =", round(log(7/6) / log(1.5), 4))
print("lambdas > 0 on (0,1) (grid):", all(f(c, t) > 0 for c in (lf, lp, lA) for t in ts[1:-1]))
