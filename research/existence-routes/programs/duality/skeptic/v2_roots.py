# v-ex2-duality session 2: independent Theta-set computation for every relation the claim quotes.
# Theta_rho = {t in [0,1] : D(t) = sum_paid p r^t - sum_targets w k^t >= 0}; reported as a union of intervals,
# endpoints by bisection to 1e-13 on a sign-change grid (D is an exponential sum; Descartes bounds its zeros).
from fractions import Fraction as F
def D(rule, t):
    T, P = rule
    return sum(p * r ** t for r, p in P.items()) - sum(w * k ** t for k, w in T.items())
def theta_set(rule, grid=4000):
    pts = [i / grid for i in range(grid + 1)]
    vals = [D(rule, t) for t in pts]
    ivs, start = [], None
    def root(a, b):
        fa = D(rule, a)
        for _ in range(80):
            mid = (a + b) / 2
            if (D(rule, mid) >= 0) == (fa >= 0): a = mid
            else: b = mid
        return (a + b) / 2
    for i, v in enumerate(vals):
        ok = v >= -1e-14
        if ok and start is None: start = pts[i] if i == 0 or vals[i-1] >= -1e-14 else root(pts[i-1], pts[i])
        if not ok and start is not None:
            ivs.append((start, root(pts[i-1], pts[i]))); start = None
    if start is not None: ivs.append((start, 1.0))
    return ivs
def inter(a, b):
    out = []
    for (x1, y1) in a:
        for (x2, y2) in b:
            lo, hi = max(x1, x2), min(y1, y2)
            if lo <= hi + 1e-15: out.append((lo, hi))
    return out
def signchanges(rule):
    T, P = rule
    coeff = {}
    for r, p in P.items(): coeff[r] = coeff.get(r, 0) + p
    for k, w in T.items(): coeff[k] = coeff.get(k, 0) - w
    seq = [coeff[b] for b in sorted(coeff) if coeff[b] != 0]
    return sum(1 for i in range(len(seq) - 1) if (seq[i] > 0) != (seq[i + 1] > 0))
rules = {
 "1+4 in 10 {3:1,2:4}<={1:10}": ({3: 1, 2: 4}, {1: 10}),
 "grouped test {3:1,2:4}<={1:6,2:2}": ({3: 1, 2: 4}, {1: 6, 2: 2}),
 "recycled {3:1,2:2}<={1:6}": ({3: 1, 2: 2}, {1: 6}),
 "grouped mixed {4:1,3:2}<={1:5,2:2}": ({4: 1, 3: 2}, {1: 5, 2: 2}),
 "pure grouped 2[3]<=3[1]+[2]": ({3: 2}, {1: 3, 2: 1}),
 "R1 2[3]<=5[1]": ({3: 2}, {1: 5}),
 "R2 2[2]<=1.05[1]+[3]": ({2: 2}, {1: 1.05, 3: 1}),
 "R2b 2[4]<=[1]+[7]": ({4: 2}, {1: 1, 7: 1}),
 "Rs 2[2]<=[1]+[3]": ({2: 2}, {1: 1, 3: 1}),
}
S = {}
for name, r in rules.items():
    S[name] = theta_set(r)
    print(f"{name}: sign changes {signchanges(r)}; Theta = {[(round(a, 6), round(b, 6)) for a, b in S[name]]}")
import math
print("log_3(5/2) =", round(math.log(2.5) / math.log(3), 6))
for other in ["R2 2[2]<=1.05[1]+[3]", "R2b 2[4]<=[1]+[7]", "Rs 2[2]<=[1]+[3]"]:
    I = inter(S["R1 2[3]<=5[1]"], S[other])
    print(f"R1 + {other}: intersection {[(round(a, 6), round(b, 6)) for a, b in I]}; max {round(max(b for a, b in I), 6)}")
# identical D for grouped test vs recycled rule (exact, as polynomials in 2^t, 3^t): coefficients
g, rr = rules["grouped test {3:1,2:4}<={1:6,2:2}"], rules["recycled {3:1,2:2}<={1:6}"]
def co(rule):
    T, P = rule; c = {}
    for r, p in P.items(): c[r] = c.get(r, 0) + p
    for k, w in T.items(): c[k] = c.get(k, 0) - w
    return {b: v for b, v in sorted(c.items()) if v != 0}
print("D coefficients grouped test:", co(g), " recycled:", co(rr), " identical:", co(g) == co(rr))
# hub what-if: paid ranks vs target ranks (single crossing?)
paid = {1, 2, 3, 4, 5, 8}; targets = {73 - r for r in (1, 2, 3, 7, 16, 48)} | {73}
print("hub what-if: max paid rank", max(paid), " min target rank", min(targets), " single crossing:", max(paid) < min(targets))
# hook / start-lift / W2 certificate profiles on Rs (paid [1]+[3] vs targets 2[2]) at block 3
for name, (p1, p2, p3) in {"hook 1:1.8:2.8": (1, 1.8, 2.8), "W2 lift 41/35:1307/630:3": (41/35, 1307/630, 3),
                            "cert A 4:6:7": (4, 6, 7), "cert B 3:5:6": (3, 5, 6), "cert C 6:7:7": (6, 7, 7),
                            "cert D 5:6:6": (5, 6, 6), "cert E 12:12:15": (12, 12, 15)}.items():
    print(f"Rs under {name}: paid {p1 + p3:.4f} vs targets {2 * p2:.4f} -> {'allowed' if p1 + p3 >= 2 * p2 - 1e-12 else 'EXCLUDED'}")
