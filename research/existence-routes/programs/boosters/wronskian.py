# e3-boost wronskian.py: numerical check of the two Wronskian closed forms in NOTE section 3 (Corollary B2).
from math import log
a, b = log(2), log(3)
lf = lambda t: 4 - 5 * 2 ** t + 2 * 3 ** t
lp = lambda t: 1 - 3 * 2 ** t + 2 * 3 ** t
lA = lambda t: -1 + 2 * 2 ** t - 3 ** t
d = lambda f, t, h=1e-6: (f(t + h) - f(t - h)) / (2 * h)
mx = 0
for i in range(1, 1000):
    t = i / 1000; x, y = 2 ** t, 3 ** t
    WA = d(lA, t) * lp(t) - lA(t) * d(lp, t); WA2 = x * y * (a * (1 - 3 ** -t) - b * (1 - 2 ** -t))
    Wf = d(lf, t) * lp(t) - lf(t) * d(lp, t); Wf2 = x * y * (7 * a * 3 ** -t - 6 * b * 2 ** -t + 4 * (b - a))
    mx = max(mx, abs(WA - WA2), abs(Wf - Wf2))
print("max |Wronskian - closed form| on a 999-point grid of (0,1):", mx)
print("h(0) =", 3 * a - 2 * b, " h(1) =", -5 * a / 3 + b, " h' = 0 at theta =", log(7 / 6) / log(1.5))
