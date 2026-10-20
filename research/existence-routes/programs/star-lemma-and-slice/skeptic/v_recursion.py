#!/usr/bin/env python3
"""v_recursion.py (v-e3-slice6): my own star recursion, written from the lemma statement only.
Local problem L(n,s): E nonsingular on AG(n,2), rank <= s on every affine hyperplane. Star lemma with Psi a d-flat,
g = 2^(n-d)-1 pieces, sigma = bound on (d+1)-flats, k = rank on Psi: 2(2^d - k) <= g (sigma - k)  (a >= 0 dropped).
So k <= floor((g sigma - 2^(d+1)) / (g-2)) for g >= 3; also k <= min(2^d, sigma). Infeasible if a bound < 0, or rank 0
on points/lines (E = 0), or Lemma 11 (s_d <= d-1: any s_d+2 points lie in a d-flat).
For U_m with W arrays: n = m-2, s = 2^(m-2) - W. Report the least excluded W per m and the open-above-bar range.
"""
import math
def bounds(n, s):
    b = {n - 1: min(s, 2 ** (n - 1))}
    for d in range(n - 2, -1, -1):
        g = 2 ** (n - d) - 1
        sig = b[d + 1]
        num = g * sig - 2 ** (d + 1)
        if num < 0: b[d] = -1; return b, 'negative at d=%d' % d
        k = min(num // (g - 2), 2 ** d, sig)
        b[d] = k
        if k == 0: return b, 'rank 0 on %d-flats' % d
        if d >= 1 and k <= d - 1: return b, 'Lemma 11 at d=%d (s_d=%d)' % (d, k)
    return b, None
def saving(m, W):
    M = 2 ** m - 1
    return math.log(m * W / M) / math.log(m)
for m in range(4, 15):
    n = m - 2; full = 2 ** n
    exW = None
    for W in range(1, full + 1):
        b, why = bounds(n, full - W)
        if why is not None and exW is None: exW = W; reason = why
        if exW is not None and why is None: print("  NON-MONOTONE at m", m, "W", W)
    barW = min([W for W in range(1, full + 1) if saving(m, W) >= 0.1] or [full + 1])
    openr = (barW, exW - 1) if barW <= exW - 1 else None
    print("m=%2d n=%2d full=%5d  excluded W >= %5d (fraction %.4f, %s)  bar W >= %d  open above bar: %s" % (
        m, n, full, exW, exW / full, reason, barW,
        ("%d..%d (savings %.3f..%.3f)" % (openr[0], openr[1], saving(m, openr[0]), saving(m, openr[1]))) if openr else "none"))
print("check m=6 W=13:", bounds(4, 3)); print("check m=6 W=12:", bounds(4, 4)); print("m=7 W=25:", bounds(5, 7), " W=26:", bounds(5, 6))
