# pf-induction: checks quoted in PROOF.md section 9 and section 5.
# (1) Lemma 1(a) on the projected block-5 and block-6 profiles (R_(r+s) <= R_r + R_s, R_r <= r, nondecreasing);
# (2) the limiting law of the projected distance at block 4 for an opposite pair from block 60 (weights on f_4..f_1).
import sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nrs_closed import N_closed
def step(a, m):
    return [0] + [sum(c * a[s] for s, c in N_closed(m, r).items() if 1 <= s <= m - 1) for r in range(1, m + 1)]
a = [0, 27456, 48576, 68320, 86592]
print("base block-4 profile ratios", [round(x / a[1], 4) for x in a[1:]],
      "subadditive:", all(a[r + s] <= a[r] + a[s] for r in range(1, 5) for s in range(1, 5) if r + s <= 4))
for m in (5, 6):
    a = step(a, m); R = [F(x, a[1]) for x in a]
    sub = all(R[r + s] <= R[r] + R[s] for r in range(1, m + 1) for s in range(1, m + 1) if r + s <= m)
    print(f"block {m}: profile {[round(float(x), 4) for x in R[1:]]}; nondecreasing {all(R[r] <= R[r+1] for r in range(1, m))}; "
          f"R_(r+s) <= R_r + R_s: {sub}; R_r <= r: {all(R[r] <= r for r in range(1, m + 1))}")
# (2) distribution of D_4 from D_60 = 60 (exact chain)
dist = {60: F(1)}
for k in range(60, 4, -1):
    P = 4 ** k - 1; new = {}
    for r, pr in dist.items():
        if r == 0: new[0] = new.get(0, 0) + pr; continue
        for s, c in N_closed(k, r).items():
            new[s] = new.get(s, 0) + pr * F(c, P)
    dist = new
print("law of the block-4 distance of a projected opposite pair (from block 60):",
      {s: round(float(p), 4) for s, p in sorted(dist.items(), reverse=True)})
S = F(1)
for k in range(5, 61): S *= F(4 ** k - 3 * 2 ** (k - 1) - 1, 4 ** k - 1)
f = {0: 0, 1: 27456, 2: 48576, 3: 68320, 4: 86592}
w = {s: float(dist.get(s, 0) / S) for s in (4, 3, 2, 1)}
print("T_60 = sum_s w_s R_s with w_s = P(D_4 = s)/S_60:", {s: round(x, 4) for s, x in w.items()},
      "->", round(sum(w[s] * f[s] / f[1] for s in w), 6))
