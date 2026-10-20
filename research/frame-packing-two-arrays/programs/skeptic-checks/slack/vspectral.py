#!/usr/bin/env python3
"""v-slack: spectral ceiling of the cut ratio on the dual polar graph (all cuts, not only point cuts), m = 2, 3, 4.
ratio(S) <= max_j (1 - P_m(j)/k_m) / (1 - theta_j/k_1), from the distance-regular structure (intersection numbers
counted directly on the enumerated Lagrangians). Eigenvalues theta_j found as the integer roots of the
characteristic polynomial of the tridiagonal intersection matrix."""
import sys
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vframes import lagrangians, dist
for m in (2, 3, 4):
    L0, Ls = lagrangians(m)
    N = len(Ls)
    d0 = [dist(m, L0, L) for L in Ls]
    # intersection numbers from L0 and one vertex at each distance
    c, a, b = {}, {}, {}
    for i in range(m + 1):
        v = next(L for L, dd in zip(Ls, d0) if dd == i)
        nb = [L for L in Ls if dist(m, v, L) == 1]
        dn = [dist(m, L0, L) for L in nb]
        c[i] = dn.count(i - 1); a[i] = dn.count(i); b[i] = dn.count(i + 1)
    k = b[0]
    # tridiagonal matrix T (rows: distance classes); eigenvalues = integer roots in [-k, k]
    def charpoly_det(x):
        # determinant of (T - xI) by recursion on the tridiagonal structure
        # T[i][i]=a_i, T[i][i+1]=b_i, T[i+1][i]=c_{i+1}
        p_prev, p = 1, (a[0] - x)
        for i in range(1, m + 1):
            p_prev, p = p, (a[i] - x) * p - b[i - 1] * c[i] * p_prev
        return p
    thetas = [x for x in range(-k, k + 1) if charpoly_det(x) == 0]
    def vpoly(x):
        v = [F(1), F(x)]
        for i in range(1, m):
            v.append(((x - a[i]) * v[i] - b[i - 1] * v[i - 1]) / F(c[i + 1]))
        return v
    km = vpoly(k)[m]
    ratios = {th: (1 - vpoly(th)[m] / km) / (1 - F(th, k)) for th in thetas if th != k}
    print('m=%d N=%d k=%d k_m=%s c=%s a=%s b=%s thetas=%s' % (m, N, k, km, c, a, b, thetas))
    print('   per-eigenspace ratio bound:', {t: str(r) for t, r in ratios.items()}, ' max =', str(max(ratios.values())),
          ' 2(1-2^-m) =', str(2 * (1 - F(1, 2 ** m))))
