"""x-pack: PROOF (exact enumeration) that no bank table of Shea's C588 residual list has a cheaper pad invoice than
P* = 724 g(1) + g(37)  (5,009 banks: 724 banks (9 x 8) with pad 1, one bank (9 x 4) with pad 37, the rest full),
for the price g_a(s) = s ((73/s)^a - 1) at a in [0.00123, 0.00125].

Facts used (all necessary conditions for ANY table that passes GX.schedK):
 (F1) per bank: n items of widths in {4,5,7,8,9}, load 73 - s with 0 <= s < 73, so 4n <= 73 - s <= 9n;
 (F2) invariant: sum over banks of w_b = 73 (9 - n_b) - 9 s_b equals 9 T - 73 I = 46368 (T = 72*5068, I = 72*616);
      banks with s = 0 have w = 73 (9 - n) <= 0 (n >= 9 by F1), so pad banks have w-sum 46368 + 73 K, K >= 0.
Let N1 = number of pad banks with (s, n) = (1, 8) (w = 64).  The others ("O") each cost >= gmin.
 * N1 >= 731: cost >= 731 g(1) > P*.
 * N1 <= N1lo: O must carry w >= 46368 - 64 N1; each O bank has w <= rho g(s) (rho = best ratio over O types),
   so cost >= N1 g(1) + (46368 - 64 N1)/rho > P* (checked below; the bound decreases in N1).
 * N1lo < N1 <= 730: O has at most kmax = floor((P* - N1 g(1)) / gmin) banks; exact DP over k <= kmax banks of the
   minimal O-cost for every exact w-sum; feasible iff 64 N1 + wsum >= 46368 and = 46368 (mod 73)."""
import sys, math
M = 73; T = 72 * 5068; I = 72 * 616; TARGET = 9 * T - 73 * I
assert TARGET == 46368

def check(a):
    def g(s): return 0.0 if s == 0 else s * ((M / s) ** a - 1)
    P = 724 * g(1) + g(37)
    O = []
    for s in range(1, M):
        for n in range(1, 19):
            if 4 * n <= M - s <= 9 * n and (s, n) != (1, 8):
                O.append((M * (9 - n) - 9 * s, g(s), s, n))
    gmin = min(c for _, c, _, _ in O)
    rho = max(w / c for w, c, _, _ in O if w > 0)
    g1 = g(1)
    # N1 <= N1lo region
    # the crude bound N1 g1 + (TARGET - 64 N1)/rho is decreasing in N1 (slope g1 - 64/rho < 0): every N1 <= N1lo
    # (the largest N1 where it still exceeds P*) is excluded
    assert g1 - 64 / rho < 0
    N1lo = 0
    while (N1lo + 1) * g1 + (TARGET - 64 * (N1lo + 1)) / rho > P * (1 + 1e-12):
        N1lo += 1
    assert 731 * g1 > P * (1 + 1e-12)          # N1 >= 731 excluded
    best = (P, 724, None)
    for N1 in range(N1lo + 1, 731):
        budget = P * (1 + 1e-12) - N1 * g1
        if budget < 0: continue
        kmax = int(budget / gmin)
        # DP: dict wsum -> min cost, over exactly k banks, k = 0..kmax
        layer = {0: 0.0}
        allbest = dict(layer)
        for k in range(1, kmax + 1):
            nxt = {}
            for ws, c in layer.items():
                for w, cc, s, n in O:
                    c2 = c + cc
                    if c2 > budget: continue
                    key = ws + w
                    if c2 < nxt.get(key, 1e18): nxt[key] = c2
            layer = nxt
            for ws, c in layer.items():
                if c < allbest.get(ws, 1e18): allbest[ws] = c
        for ws, c in allbest.items():
            tot = 64 * N1 + ws
            if tot >= TARGET and (tot - TARGET) % M == 0:
                cost = N1 * g1 + c
                if cost < best[0] * (1 - 1e-12):
                    best = (cost, N1, ws)
    return P, best, N1lo, gmin, rho

for a in (0.00123, 0.0012383, 0.00125):
    P, best, N1lo, gmin, rho = check(a)
    print('a=%.7f  P*/a=%.6f  N1 <= %d excluded by the ratio bound (rho=%.4f, gmin/a=%.4f);  best found /a=%.6f at N1=%d  -> %s'
          % (a, P / a, N1lo, rho / (1 / a), gmin / a, best[0] / a, best[1], 'P* OPTIMAL' if best[2] is None else 'IMPROVED %s' % (best,)))
