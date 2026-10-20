#!/usr/bin/env python3
"""t_sph_exact.py: exact check (Fractions) for m = 2..14 that theta_j = 2(2^{m-j}-1) - (2^j - 1) are the eigenvalues of
the dual polar graph's intersection matrix (b_i = 2^{i+1}(2^{m-i}-1), c_i = 2^i - 1), and that the ratios used in
Theorem E2(ii) are c_1/mu_1 = 2 - 2^{1-m} and c_m/mu_m = (2/3)(1 - (-1/2)^m), so c and mu are never proportional."""
from fractions import Fraction as Fr
ok = True
for m in range(2, 15):
    b = [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)]; c = [2 ** i - 1 for i in range(m + 1)]
    k = b[0]; a = [k - b[i] - c[i] for i in range(m + 1)]
    om = {}
    for j in range(m + 1):
        th = 2 * (2 ** (m - j) - 1) - (2 ** j - 1)
        w = [Fr(1), Fr(th, k)]
        for i in range(1, m): w.append((th * w[i] - c[i] * w[i - 1] - a[i] * w[i]) / b[i])
        last = c[m] * w[m - 1] + a[m] * w[m] - th * w[m]
        ok &= (last == 0); om[j] = w
    r1 = (om[1][m] - 1) / (om[1][1] - 1); rm = (om[m][m] - 1) / (om[m][1] - 1)
    ok &= (r1 == 2 - Fr(2, 2 ** m)) and (rm == Fr(2, 3) * (1 - Fr(-1, 2) ** m)) and (r1 != rm)
    ok &= len({2 * (2 ** (m - j) - 1) - (2 ** j - 1) for j in range(m + 1)}) == m + 1
    print("m=%2d: eigenvalue check %s; c_1/mu_1 = %s, c_m/mu_m = %s" % (m, last == 0, r1, rm))
print("ALL EXACT CHECKS PASS" if ok else "FAILED")
