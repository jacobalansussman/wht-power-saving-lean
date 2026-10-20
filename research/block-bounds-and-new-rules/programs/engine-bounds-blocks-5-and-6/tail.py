# pf-induction, 13 Oct 2026.  The ALL-m statements of PROOF.md (Theorems 4 and 5), exact arithmetic.
#  (1) exact projection recursion from the block-4 hook (no triangle improvement) up to block M = 60;
#  (2) checks that the block-M profile R_r = a_r / a_1 is nondecreasing in r (needed for the coupling);
#  (3) the dominating deficit chain: up 2^-(j+1), down 1 - 2^-j (j >= 1); checks p_up(j,k) <= 2^-(j+1) and
#      p_down(j,k) >= 1 - 2^-j exactly for all 1 <= k <= 80, 0 <= j < k;
#  (4) its stationary law pi*(j) ~ 1 / prod_{i<=j} (2^i - 1) and the all-m lower bound
#      c* = sum_j pi*(j) R^(M)_(M-j)  (normalizer bounded ABOVE, so c* is a certified lower bound);
#  (5) the no-growth upper bound for ANY base profile: T_m <= max_s R_s / S_inf(m0), S = prod_k A_1(k)/(4^k - 1).
# Usage: python3 -I tail.py
import sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nrs_closed import N_closed

def step(a, m):
    b = [0] * (m + 1)
    for r in range(1, m + 1):
        b[r] = sum(c * a[s] for s, c in N_closed(m, r).items() if 1 <= s <= m - 1)
    return b

def project(base, M):
    a = [0] + list(base); m0 = len(base); tops = {}
    for m in range(m0 + 1, M + 1):
        a = step(a, m); tops[m] = F(a[m], a[1])
    return a, tops

M = 60
bases = {"hook m=4 k=8 (27456, 48576, 68320, 86592)": (27456, 48576, 68320, 86592),
         "hook m=3 k=4 (80, 144, 224)": (80, 144, 224),
         "perfect m=2 (1, 2)": (1, 2)}

# (3) dominating chain inequalities, exact
viol = 0
for k in range(1, 81):
    P = 4 ** k - 1
    for j in range(0, k):
        r = k - j
        Nr = N_closed(k, r)
        pup = F(Nr.get(r - 2, 0), P); pdown = F(Nr.get(r, 0), P)
        if pup > F(1, 2 ** (j + 1)): viol += 1
        if j >= 1 and pdown < 1 - F(1, 2 ** j): viol += 1
print(f"(3) dominating chain: p_up(j,k) <= 2^-(j+1) and p_down(j,k) >= 1 - 2^-j for all k <= 80, j < k: violations {viol}")
# monotonicity condition of the dominating chain: p*_up(j) + p*_down(j+1) <= 1
mono = all(F(1, 2 ** (j + 1)) + (1 - F(1, 2 ** (j + 1))) <= 1 for j in range(0, 200))
print(f"    dominating chain is stochastically monotone (p_up(j) + p_down(j+1) <= 1): {mono}")

# (4) stationary law, exact on j <= J, tail of the normalizer bounded above
J = 80
w = [F(1)]
for j in range(1, J + 1): w.append(w[-1] / (2 ** j - 1))
Zlo = sum(w); tail_up = 2 * w[J]          # w_{j+1}/w_j <= 1/3 for j >= 1, so sum_{j>J} w_j <= w_J/2 < 2 w_J
Zup = Zlo + tail_up
print(f"    pi*(0..4) = {[round(float(x / Zlo), 6) for x in w[:5]]}, normalizer in [{float(Zlo):.12f}, {float(Zup):.12f}]")
for name, base in bases.items():
    a, tops = project(base, M)
    R = [None] + [F(a[r], a[1]) for r in range(1, M + 1)]
    nondecr = all(R[r] <= R[r + 1] for r in range(1, M))
    cstar = sum(w[j] * R[M - j] for j in range(0, min(J, M - 1) + 1)) / Zup    # R_(M-j) >= 0; dropped terms >= 0
    mn = min(tops.values()); argmn = min(tops, key=lambda m: tops[m])
    print(f"== {name}")
    print(f"   exact tops: " + ", ".join(f"m={m}: {float(tops[m]):.6f}" for m in sorted(tops) if m <= 8 or m in (10, 20, 40, 60)))
    print(f"   min over m = {min(tops)}..{M} of the exact top ratio: {float(mn):.9f} (at m={argmn})")
    print(f"   (2) block-{M} profile nondecreasing in r: {nondecr}")
    print(f"   (4) certified bound for EVERY m >= {M}: T_m >= c* = {float(cstar):.9f}")
    # (5) no-growth upper bound
    m0 = len(base)
    S = F(1)
    for k in range(m0 + 1, 400): S *= F(4 ** k - 3 * 2 ** (k - 1) - 1, 4 ** k - 1)
    # remaining factors k >= 400 change S by less than 1e-100; bound: prod_{k>=400}(1 - 2^-k) >= 1 - 2^-399
    Slow = S * (1 - F(1, 2 ** 399))
    Rmax = max(F(x, base[0]) for x in base)
    print(f"   (5) no-growth: S_inf(m0={m0}) >= {float(Slow):.9f}; every block-m projection of this base has T_m <= {float(Rmax / Slow):.6f}")
    # limit of the projection, for reference
    print(f"   limit of the exact tops (m = {M}): {float(tops[M]):.12f}")
# (5') the generic ceiling of projection from block m0 with ANY family (R_max <= m0):
for m0 in range(2, 9):
    S = F(1)
    for k in range(m0 + 1, 400): S *= F(4 ** k - 3 * 2 ** (k - 1) - 1, 4 ** k - 1)
    print(f"   projection from block {m0}, ANY base family: T_m <= {m0} / S_inf = {float(m0 / (S * (1 - F(1, 2 ** 399)))):.6f} for every m")

# (6) the headline constant, compared exactly
base = bases["hook m=4 k=8 (27456, 48576, 68320, 86592)"]
a, tops = project(base, M)
R = [F(0)] + [F(a[r], a[1]) for r in range(1, M + 1)]
cstar = sum(w[j] * R[M - j] for j in range(0, min(J, M - 1) + 1)) / Zup
mn = min(tops.values())
print(f"(6) block-4 hook: c* = {cstar.numerator * 10**15 // cstar.denominator / 10**15:.15f}, "
      f"min exact top (m=5..{M}) = {mn.numerator * 10**15 // mn.denominator / 10**15:.15f}; "
      f"all-m bound = min of the two >= 2.9447: {min(cstar, mn) >= F(29447, 10000)}")
print(f"    exact: T_5 = {tops[5]}, T_6 = {tops[6]}, T_7 = {tops[7]}")
