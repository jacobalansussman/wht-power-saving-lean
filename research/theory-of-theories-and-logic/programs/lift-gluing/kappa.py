# Independent check of kappa_m: least weight fraction of a set of levels i in 1..m-2 meeting every adjacent pair
# (i, i+1), i = 0..m-2, where levels 0 and m-1 may not be used (phi_0 = phi_{m-1} = 0). Weight of level i:
# e_{m,i} = [m,i]_2 (2^(m-i) - 1). Exact rationals; DP over the path.
from fractions import Fraction as Fr
def gb(m, k):
    num = den = 1
    for i in range(k): num *= (2 ** (m - i) - 1); den *= (2 ** (i + 1) - 1)
    return num // den
def kappa(m):
    w = [gb(m, i) * (2 ** (m - i) - 1) for i in range(m)]; tot = sum(w)
    INF = None
    # dp[in] = least weight of a valid cover of levels 0..i with level i in (1) / out (0)
    dp = {0: 0, 1: INF}       # level 0 is out
    for i in range(1, m):
        allowed = 1 <= i <= m - 2
        new = {}
        new[0] = dp[1]        # i out -> i-1 must be in
        new[1] = (w[i] + min(v for v in dp.values() if v is not None)) if allowed else INF
        dp = new
    best = min(v for v in dp.values() if v is not None)
    return Fr(best, tot), tot
# also the linear-program bound of Lemma L at nu = 1 with phi >= 0 continuous (should equal integral VC / 3)
for m in list(range(3, 13)) + [20, 21, 40, 41, 80, 81, 160, 161]:
    k, tot = kappa(m)
    print(m, "e_m =", tot if m < 13 else "...", "kappa = %.6f" % float(k), " (2/3)kappa = %.6f" % (2 * float(k) / 3))
for m in (3, 4, 5, 6, 8):
    w = [gb(m, i) * (2 ** (m - i) - 1) for i in range(m)]
    print("m", m, "levels", w, "e_m", sum(w))
print("thresholds m^0.1-1:", [(m, round(m ** 0.1 - 1, 4)) for m in (3, 4, 5, 6, 17, 18)])
print("canonical m=8 per edge 2*3927388/e_8 =", 2 * 3927388 / sum(gb(8, i) * (2 ** (8 - i) - 1) for i in range(8)))
print("decomposed per edge: 2*1/35, 36/240, 580/2077, 2*5384/e6:", 2/35, 36/240, 580/2077, 2*5384/sum(gb(6,i)*(2**(6-i)-1) for i in range(6)))
