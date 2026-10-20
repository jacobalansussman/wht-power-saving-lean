#!/usr/bin/env python3
"""ex2-reform calibration: known exponents re-derived as theta-price (Collatz-Wielandt) thresholds,
and the 'most favourable transfer' of OpenAI's reflection contradiction to theta-prices.

Part A. One-type rule: W tasks of sizes m_j (W_j of each) from calls of sizes r_t.
        theta-price feasible  <=>  sum_j W_j m_j^theta <= sum_t r_t^theta   (v cancels).
        The threshold theta* solves equality; saving = 1 - theta*.  This is Theorem C's z,
        ex-laser's tau-transfer, and OpenAI's amplification exponent (g rho^alpha = 1).
Part B. Two types with cross rules: nonlinear Perron root of T(v)_A = min_rules sum rho^theta v_B.
        Re-derives universality (iso 2b): a strongly connected price system has one threshold.
Part C. OpenAI section 8 inequalities with the tensor law replaced by the power law (theta-world):
        lower: pi(R) >= (2d)^theta/2 - 4 ;  upper: pi(R) <= beta(theta)*(d/2+1) + 2,
        beta(theta) = 3*2^(-theta) (their dense-pair bound uses p(H(x)H)=4p(H)).
        Contradiction at d  <=>  lower > upper.  At theta = 1 this must give d > 30 (paper).
"""
import math

def root(f, lo, hi, it=200):
    flo = f(lo)
    for _ in range(it):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2

print("== Part A: one-type thresholds (calibration against recorded values) ==")
def one_type(targets, calls):
    # targets: list of (count, size); calls: list of sizes. Returns theta* with sum W m^th = sum r^th
    f = lambda th: sum(W * m ** th for W, m in targets) - sum(r ** th for r in calls)
    return root(f, 1e-9, 1.0)

cases = [
    ("2 in 3 at block 2 (rank-1 moves)", [(2, 2)], [1] * 3, "0.415 (ex-note)"),
    ("3 in 5 at block 2", [(3, 2)], [1] * 5, "0.263 (ex-note)"),
    ("2 arrays block 3, ranks 1,1,1,1,1", [(2, 3)], [1] * 5, "x = 0.1660 (w-engine 4.2)"),
    ("2 arrays block 3, ranks 2,1,1,1", [(2, 3)], [2, 1, 1, 1], "x = 0.2187 (w-engine 4.2)"),
    ("2 arrays block 3, ranks 2,2,1", [(2, 3)], [2, 2, 1], "x = 0.3272 (w-engine 4.2)"),
    ("3 arrays block 3, eight rank 1", [(3, 3)], [1] * 8, "x = 0.1072 (w-engine 4.2)"),
    ("3 arrays block 3, ranks 2,2,2,2", [(3, 3)], [2, 2, 2, 2], "x = 0.2905 (w-engine 4.2)"),
    ("mixed 1 + 4 in 10 at block 3 (tau-transfer)", [(1, 3), (4, 2)], [1] * 10, "0.119 (ex-laser Thm 1)"),
    ("mixed 7 + 38 in 88 at block 3", [(7, 3), (38, 2)], [1] * 88, "0.125 (ex-laser Thm 1)"),
    ("5 in 14 at block 3", [(5, 3)], [1] * 14, "0.063 (ex-note)"),
    ("2 in 7 at block 4", [(2, 4)], [1] * 7, "0.096 (meta2)"),
    ("3 in 10 at block 4", [(3, 4)], [1] * 10, "0.132 (meta2)"),
    ("grouped 2 arrays block 4, ranks 2,2,2,1", [(2, 4)], [2, 2, 2, 1], "0.1675 (MAP 0.5)"),
    ("plain W in Wm at block 5 (must give saving 0)", [(4, 5)], [1] * 20, "0 (plain FFT)"),
]
for name, T, C, rec in cases:
    th = one_type(T, C)
    print(f"  {name:48s} theta* = {th:.4f}  saving = {1-th:.4f}   recorded: {rec}")

# OpenAI amplification: g calls, rho = 1/(b q^(b-1)); alpha = ln g / ln(b q^(b-1)).  Same as one-type rule
# with one target of 'size' b q^(b-1) and g calls of size 1.
for (q, b, g) in [(2, 2, 3), (2, 3, 11), (2, 4, 31)]:
    S = b * q ** (b - 1)
    th = one_type([(1, S)], [1] * g)
    print(f"  OpenAI-type win q={q} b={b} g={g} (plain {S}): alpha = {th:.4f} = ln g/ln(b q^(b-1)) = {math.log(g)/math.log(S):.4f}")

print()
print("== Part B: two types, cross rules (universality as strong connectivity) ==")
# Types H and U. Rules (normalized cost per coordinate ~ v * n^theta):
#   H unit: 2 in 3 at block 2 (a stand-in unit; any saving unit for H would do): 2*2^th*vH <= 3*vH
#   U from 3 calls of H at the same size (shrink 1): vU <= 3 vH
#   H from 2 calls of U at the same size: vH <= 2 vU
# Nonlinear Perron: T(v)_H = min(3*2^-th*vH/2 ... ) -- with the H-unit the H coordinate contracts
# whenever 3/(2*2^th) < 1.  U's threshold equals H's because U -> H -> (unit) -> H ... 
def perron_two(th, iters=3000):
    # normalized price operator: for each type, min over its rules of (sum over calls of size^th * v) / (target size^th)
    # H rules: plain block 2 (2 calls of size 1 for a size-2 target): 2*2^-th*vH ; unit 2-in-3: 1.5*2^-th*vH ;
    #          H from 2 calls of U at the same size: 2*vU.
    # U rules: plain block 2: 2*2^-th*vU ; U from 3 calls of H at the same size: 3*vH.
    vH, vU = 1.0, 1.0
    for _ in range(iters):
        nH = min(2 * 2 ** -th * vH, 1.5 * 2 ** -th * vH, 2 * vU)
        nU = min(2 * 2 ** -th * vU, 3 * vH)
        s = max(nH, nU)
        lam = s / max(vH, vU)
        vH, vU = nH / s, nU / s
    return lam, vH, vU
for th in [0.5, 0.58, 0.585, 0.59, 0.7]:
    lam, vH, vU = perron_two(th)
    print(f"  theta={th:.3f}: cone spectral radius of the price operator = {lam:.4f}; eigvec (vH,vU)=({vH:.3f},{vU:.3f})  [<1: theta-saving for both types]")
print("  threshold for H alone: 3 = 2*2^th -> th =", f"{math.log(1.5)/math.log(2):.4f}", "(saving 0.415); U gets the same threshold")

print()
print("== Part C: OpenAI's reflection contradiction, most favourable transfer to theta-prices ==")
def gap(d, th, beta=None, low_b=4.0, up_b=2.0):
    if beta is None:
        beta = 3.0 * 2 ** (-th)
    lower = (2 * d) ** th / 2 - low_b
    upper = beta * (d / 2 + 1) + up_b
    return lower - upper
# calibration at theta = 1
dmin = next(d for d in range(1, 10 ** 6) if gap(d, 1.0) > 0)
print(f"  theta = 1: first d with a contradiction = {dmin}  (paper: incompatible when d > 30)")
def best(th, beta=None, dmax=10 ** 7):
    bd, bg = None, -1e18
    d = 1
    while d <= dmax:
        g = gap(d, th, beta)
        if g > bg:
            bd, bg = d, g
        d = d + 1 if d < 2000 else int(d * 1.01) + 1
    return bd, bg
f = lambda th: best(th)[1]
thc = root(lambda th: f(th), 0.80, 1.0, it=60)
bd, bg = best(thc + 1e-6)
print(f"  critical theta (largest theta with NO contradiction at any d): theta_c = {thc:.4f}  -> saving at most {1-thc:.4f}; best d there ~ {bd}")
for th in [0.99, 0.97, 0.96, 0.95, 0.9]:
    d, g = best(th)
    print(f"  theta = {th}: max over d of (lower - upper) = {g:+.3f} at d = {d}")
# what dense-pair constant would refute theta = 0.9?
def need_beta(th):
    return root(lambda b: best(th, beta=b)[1], 0.01, 3.0, it=60)
for th in [0.9, 0.95]:
    b = need_beta(th)
    d, g = best(th, beta=b * 0.999)
    print(f"  to refute theta = {th}: dense-pair price must be < {b:.4f} (OpenAI: 3*2^-theta = {3*2**-th:.4f}); best d ~ {d}")
# with zero boundary costs, a constant slope ratio never refutes theta<1 at large d:
for th in [0.9]:
    d, g = best(th, beta=None)
    print(f"  (ratio check) at theta=0.9 the cheap/tensor slope ratio needed at d is (2d)^-0.1: d=150 -> {(300)**-0.1:.3f}, d=10^4 -> {(2e4)**-0.1:.3f}; OpenAI's is 3/4 (times 2^(1-theta))")
