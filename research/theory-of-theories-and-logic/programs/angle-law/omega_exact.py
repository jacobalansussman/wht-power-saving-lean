# Skeptic's own exact computation of the spherical functions of the dual polar graph DSp(2m,2).
# Independent of a2-vi-stability's code. Exact rationals.
from fractions import Fraction as F
import sys
from math import comb

def params(m):
    k = 2*(2**m - 1)
    b = [2**(j+1)*(2**(m-j)-1) for j in range(m+1)]
    c = [2**j - 1 for j in range(m+1)]
    a = [k - b[j] - c[j] for j in range(m+1)]
    return k, a, b, c

def omegas(m, theta):
    k, a, b, c = params(m)
    w = [F(1), F(theta, k)]
    for j in range(1, m):
        # c_j w(j-1) + a_j w(j) + b_j w(j+1) = theta w(j)
        w.append((theta*w[j] - a[j]*w[j] - c[j]*w[j-1]) / b[j])
    # closing condition at j = m: c_m w(m-1) + a_m w(m) = theta w(m)
    close = c[m]*w[m-1] + a[m]*w[m] - theta*w[m]
    return w[:m+1], close

def kd(m):
    k, a, b, c = params(m)
    ks = [1]
    for j in range(1, m+1):
        ks.append(ks[-1]*b[j-1]//c[j])
    return ks

M = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mingap = None
for m in range(1, M+1):
    k, a, b, c = params(m)
    ks = kd(m)
    N = sum(ks)
    Nexp = 1
    for i in range(1, m+1): Nexp *= (2**i + 1)
    assert N == Nexp, (m, N, Nexp)
    thetas = [2*(2**(m-i)-1) - (2**i - 1) for i in range(m+1)]
    W = []
    mult = []
    for i, th in enumerate(thetas):
        w, close = omegas(m, th)
        assert close == 0, ("not an eigenvalue", m, i)
        W.append(w)
        mi = F(N) / sum(ks[j]*w[j]**2 for j in range(m+1))
        assert mi.denominator == 1, ("non-integral multiplicity", m, i, mi)
        mult.append(int(mi))
    assert sum(mult) == N
    # Lemma 2 closed form
    y = 2**m
    for j in range(m+1):
        assert W[1][j] == (F(y+1, 2**j) - 2) / (y-1), ("lemma2", m, j)
    # omega_i(m) closed form claimed
    for i in range(m+1):
        assert W[i][m] == F((-1)**i, 2**(i*(m+1-i))), ("omega_i(m)", m, i, W[i][m])
    # T_m(j)
    for j in range(1, m+1):
        ratios = [(1 - W[i][j]) / (1 - W[i][1]) for i in range(1, m+1)]
        T = max(ratios)
        assert T == 2 - F(2, 2**j), ("T", m, j, T)
        if j >= 2:
            others = sorted(ratios, reverse=True)
            assert ratios[0] == T and ratios.count(T) == 1, ("argmax", m, j)
            gap = others[0] - others[1]
            if mingap is None or gap < mingap[0]:
                mingap = (gap, m, j, ratios.index(others[1]) + 1)
        # proved bound <8/3 check
        assert T < F(8, 3)
    # NN^T for frame-point incidence: common points at distance d = 2^(m-d) - 1
    lam = [sum((2**(m-d)-1)*ks[d]*W[i][d] for d in range(m+1)) for i in range(m+1)]
    nz = [i for i in range(m+1) if lam[i] != 0]
    npts = 2**(2*m) - 1
    if m <= 6 or m in (10, 20, 30, 40):
        print(f"m={m:2d} N={N} mult={mult if m<=6 else '...'} NN^T nonzero on E_i for i in {nz}; dims {[mult[i] for i in nz]} sum {sum(mult[i] for i in nz)} vs points {npts}")
    assert nz == [0, 1], ("incidence image", m, nz)
print("ALL CHECKS PASSED for m = 1 ..", M)
print("min gap (j >= 2) between top ratio (i=1) and runner-up:", float(mingap[0]), "=", mingap[0] if mingap[0].denominator < 10**12 else "", "at (m, j) =", mingap[1:3], "runner-up i =", mingap[3])
