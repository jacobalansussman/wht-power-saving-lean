# pf-induction, 13 Oct 2026.  kappa_m = the largest ratio  P_m(cut) / P_1(cut)  over all cuts (and all
# partitions) of the dual polar graph DSp(2m,2) (PROOF.md, Theorem 6).  Spectral bound, exact rationals:
#   kappa_m <= max over eigenvalues theta != k of (1 - u_m(theta)) / (1 - u_1(theta)),
# u_i(theta) the normalised cosine sequence from the intersection array b_i = 2^(i+1) (2^(m-i) - 1), c_i = 2^i - 1.
# Eigenvalues found exactly as the integer roots of the tridiagonal intersection matrix (all m+1 of them must be
# found, else the script stops).  Also checks the sphere sizes against engine_m.py (m = 3, 4, 5) and, at m = 2, 3,
# a brute-force Lagrangian list with a point cut.
# Usage: python3 -I cutceil.py
from fractions import Fraction as F
import itertools

def arrays(m):
    b = [2 ** (i + 1) * (2 ** (m - i) - 1) for i in range(m + 1)]
    c = [2 ** i - 1 for i in range(m + 1)]
    k = b[0]; a = [k - b[i] - c[i] for i in range(m + 1)]
    return b, c, a, k

def spheres(m):
    b, c, a, k = arrays(m); s = [1]
    for i in range(1, m + 1): s.append(s[-1] * b[i - 1] // c[i])
    return s

def cosines(m, th):
    b, c, a, k = arrays(m)
    u = [F(1), F(th, k)]
    for i in range(1, m):
        u.append((th * u[i] - c[i] * u[i - 1] - a[i] * u[i]) / b[i])
    ok = (c[m] * u[m - 1] + a[m] * u[m] == th * u[m])
    return u, ok

def eigen(m):
    b, c, a, k = arrays(m)
    ev = [th for th in range(-k, k + 1) if cosines(m, th)[1]]
    assert len(ev) == m + 1, (m, ev)
    return ev

assert spheres(3) == [1, 14, 56, 64] and spheres(4) == [1, 30, 280, 960, 1024] and spheres(5) == [1, 62, 1240, 9920, 31744, 32768]
print("sphere sizes agree with engine_m.py (m = 3, 4, 5)")

def brute(m):
    n2 = 2 * m; M = (1 << m) - 1
    def form(x, y): return (bin((x & M) & (y >> m)).count('1') + bin((y & M) & (x >> m)).count('1')) & 1
    def span(B):
        out = {0}
        for v in B:
            if v not in out: out |= {x ^ v for x in out}
        return frozenset(out)
    lags = set()
    for B in itertools.combinations(range(1, 1 << n2), m):
        if all(form(x, y) == 0 for x in B for y in B):
            S = span(B)
            if len(S) == 2 ** m: lags.add(S)
    return list(lags)

bounds = {}
for m in range(2, 13):
    ev = eigen(m); b, c, a, k = arrays(m)
    best = None
    for th in ev:
        if th == k: continue
        u, _ = cosines(m, th)
        val = (1 - u[m]) / (1 - u[1])
        if best is None or val > best[0]: best = (val, th)
    bounds[m] = best[0]
    print(f"m={m}: eigenvalues {ev}; kappa_m <= {best[0]} = {float(best[0]):.6f} (eigenvalue {best[1]})")

for m in (2, 3):
    lags = brute(m)
    v = 1
    f = [1 if v in L else 0 for L in lags]
    tot = {}
    for i, L in enumerate(lags):
        for j, Mx in enumerate(lags):
            d = m - (len(L & Mx).bit_length() - 1)
            s, cnt = tot.get(d, (0, 0)); tot[d] = (s + (f[i] - f[j]) ** 2, cnt + 1)
    ratio = F(tot[m][0], tot[m][1]) / F(tot[1][0], tot[1][1])
    print(f"brute m={m}: {len(lags)} Lagrangians; point cut {{U : v in U}}: P_m/P_1 = {ratio} = {float(ratio):.6f} "
          f"<= spectral {float(bounds[m]):.6f}: {ratio <= bounds[m]}")
